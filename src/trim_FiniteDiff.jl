include("simulate.jl")
using FiniteDifferences
using Ipopt
using GenericLinearAlgebra
fdm = central_fdm(5,1);

## Ipopt directly
function trim_FinDiff(j::Vector{Joint},GravityInInertial::Vector{Float64};
  ix::Vector{Integer}=zeros(Integer,20*(nB-1)+1),
  iy::Vector{Integer}=zeros(Integer,7*(nB-1)))
  nB = length(j) + 1; # number of bodies in the joint tree
  # ix: 0 indicates that the corresponding variable is free during optimization

  # iy: 0 indicates that the corresponding second derivative term needs to be zero
  # (specify indices for iy starting from second body)
  # iy = [1;4;5] means that the ẍ. β̈1 and β̈2 are free for j[1].RB2

  # Find which variables are free for optimization
  freeVarInd = findall(x->x==0,ix); # indices of free variables

  # find which acceleration-level terms are desired to be 0/constrained
  freeAccInd = findall(x->x==1,iy) # indices of free accelerations

  # Initial guesses
  # Initial guess is provided initial location
  x0 = []; # Initial guess
  for k = 1:length(j)
    append!(x0,j[k].RB2.x)
  end
  extFListCopy = externalFTotal(0.0,j);
  u0 = genU(extFListCopy); # Initial guess for control input

  n_xu = length(freeVarInd)-1; # Number of optimization variables
  lb = [-Inf*ones(n_xu);0.0]; ub = [Inf*ones(n_xu);Inf] # Bounds on optimization variables
  # lb[4] = 0.0;
  # lb[14+4] = 0.0;

  m = 1; # Number of nonlinear constraints (xdot = f(x,u)) for second order variables
  lb_g = [-Inf]; ub_g = [0.0];

  # extract states and control from optimization variables
  function getXU(z::AbstractArray{T}) where T<:Real
    z0 = [x0;u0[:,2:end][:];0.0];
    z0[freeVarInd] .= 0;

    Mz = zeros(length(z0), length(freeVarInd))
    for (i,ind) in enumerate(freeVarInd)
      Mz[ind,i] = 1;
    end

    zFull = z0 + Mz*z;

    X = zFull[1:14*(nB-1)];
    U = zFull[14*(nB-1)+1:end-1]
    return X,U
  end

  ## Constraint functions
  function qddFn(xuVar)
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    xdot = fxdot(x,u,j,g);
    y = Vector(undef,7*nB); # qdd
    for i=1:nB
      y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
    end
    deleteat!(y,1:7); # inertial frame has no acceleration
    deleteat!(y,freeAccInd); # free acceleration variables
    return y
  end

  function quatNormConstr(xuVar)
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    q_constr = Vector(undef, nB)
    for i=1:nB
        q_constr[i] = norm(x[14*(i-1)+4:14*(i-1)+7])^2 - 1;
    end
    return q_constr
  end

  function jointLocConstr(xuVar)
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    jLoc_constr = Matrix(undef,3,length(j))
    for i=1:length(j)
        b1_id = j[i].RB1.bodyID; b2_id = j[i].RB2.bodyID;
        x1 = x[(b1_id-1)*14+1:b1_id*14]
        x2 = x[(b2_id-1)*14+1:b2_id*14]
        j_b1 = x1[1:3] + transpose(quat2dcm(x1[4:7]))*j[i].pos1;
        j_b2 = x2[1:3] + transpose(quat2dcm(x2[4:7]))*j[i].pos2;
        jLoc_constr[:,i] = j_b1 - j_b2;
    end
    return jLoc_constr
  end

  function _nonlinear_constraints(z::AbstractArray{T}) where T<:Real
    xVar, uVar = getXU(z)
    x = [Vector{T}(j[1].RB1.x); xVar]
    u = [zeros(T,6) reshape(uVar,(6, nB-1))]
    xuVar = [x;u[:]]; # Vector of x,u

    # Dynamics constraint
    y = qddFn(xuVar);
    # Quaternion norm constraint
    q_constr = quatNormConstr(xuVar);
    # JointLocation constraint  (translation constraint)
    jLoc_constr = jointLocConstr(xuVar);

    constr = [y;q_constr;jLoc_constr[:]];

    J = sum(constr.*constr) - z[end];
    return J
  end

  function eval_f(z)
    J = z[end]
    return J
  end

  function eval_g(z, g)
    g[:] .= _nonlinear_constraints(z)
  end

  function eval_grad_f(z, grad_f)
    grad_f[:] = FiniteDifferences.grad(fdm, eval_f, z)[1];
  end

  function eval_jac_g(z, mode, rows, cols, values)
    nvar = length(z);
    ng = 1;

    if mode == :Structure
        index = 1;
        for i=1:ng
            for j=1:nvar
                rows[index] = i;
                cols[index] = j;
                index += 1;
            end
        end
    else
        values[:] = FiniteDifferences.grad(fdm,_nonlinear_constraints,z)[1];
    end

    # return jac_g
  end

  function grad_f_fn(z)
    grad_f = FiniteDifferences.grad(fdm, eval_f, z)[1];
    return grad_f;
  end

  function grad_nonlinear_constraints_fn(z)
    grad_nc = FiniteDifferences.grad(fdm, _nonlinear_constraints, z)[1];
    return grad_nc;
  end

  function eval_h(z, mode, rows, cols, obj_factor, lambda, values) # Hessian of cost.
    nvar = length(z);
    if mode == :Structure
      # Symmetric matrix, fill the lower left triangle only
      idx = 1
      for row = 1:nvar
        for col = 1:row
          rows[idx] = row
          cols[idx] = col
          idx += 1
        end
      end
    else
      H = obj_factor*FiniteDifferences.jacobian(fdm, grad_f_fn, z)[1];

      for i=1:length(lambda)
        H += lambda[i]*FiniteDifferences.jacobian(fdm, grad_nonlinear_constraints_fn, z)[1];
      end

      # Copy them into return variable -- only lower triangle.
      idx = 1;
      for row = 1:nvar
        for col = 1:row
          values[idx] = H[row,col];
          idx += 1;
        end
      end
    end
  end

  # Initial guess from freeVarInd
  n = length(freeVarInd);
  guessFull = [x0;u0[:,2:end][:];0.0];
  inGuess = guessFull[freeVarInd];

  nele_jac = Integer(m*n);
  nele_hess = Integer(n*(n+1)/2);
  prob = createProblem(n, lb, ub, m, lb_g, ub_g, nele_jac, nele_hess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);
  prob.x = inGuess; # Initial guess.
  status = solveProblem(prob);

  trim_x = prob.x[1:length(x0)];
  gam = prob.x[end];
  trim_u = prob.x[length(x0)+1:end-1];
  return trim_x, trim_u, gam
end
## testing the trim function
# x0Orig, u0Orig = getXU_0(j);
# nB = length(j) + 1;
# ix= zeros(Integer,20*(nB-1)+1);
# # WITH U=0
# ix[14*(nB-1)+1:20*(nB-1)] .= 1;
# trim_U = u0Orig;
# iy = zeros(Integer,7*(nB-1));
#
# out = trim_FinDiff(j,g,ix=ix, iy=iy);
# trim_x, trim_u, gam = out;
# println("trim_x = ", trim_x)
# trim_X = [j[1].RB1.x;trim_x]
# # trim_U = [zeros(6) reshape(trim_u,(6,nB-1))];
# # println("trim_U = ", trim_U[:,2])
# # trim_U = zeros(6,nB);
# println("qconstr = ", norm(trim_x[4:7])-1)
# println("ẍ =", norm(fxdot(trim_X,trim_U, j, g)))
