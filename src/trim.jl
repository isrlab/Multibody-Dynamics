include("simulateDiff.jl")

using Ipopt
using GenericLinearAlgebra
clearconsole();
nB = length(j) + 1;
## JumP
# using JuMP
# model = Model();
# set_optimizer(model, Ipopt.Optimizer)
# @variable(model, x[1:14*nB])
# @variable(model, u[1:6,1:nB])
# # @NLexpression(model, )
# obj(u...) = sum(u.^2)# for i in 1:length(u)
# register(model, :obj, length(u), obj, autodiff = true)
#
# function fxdot_constr(x..., u...)
#   xdot = fxdot(x,u,j,GravityInInertial)
#   y = zeros(7*nB);
#   y[1:7] = xdot[8:14];
#   for i=2:nB
#     y[7*(i-1)+1:7*(i)] = xdot[7*(i+1)+1:14*i]
#   end
#
#   return y
# end
# register(model, :fxdot, length(x) )
# @NLobjective(model, Min, obj(u...))

## Ipopt directly
function trimDiff(j::Vector{Joint},GravityInInertial::Vector{Float64};
  ix::Vector{Integer}=zeros(Integer,20*(nB-1)+1),
  iy::Vector{Integer}=zeros(Integer,7*(nB-1)))
  # ix: 0 indicates that the corresponding variable is free during optimization

  # iy: 0 indicates that the corresponding second derivative term needs to be zero
  # (specify indices for iy starting from second body)
  # iy = [1;4;5] means that the ẍ. β̈1 and β̈2 are free for j[1].RB2

  # Find which variables are free for optimization
  freeVarInd = findall(x->x==0,ix); # indices of free variables

  # find which acceleration-level terms are desired to be 0
  freeAccInd = findall(x->x==1,iy) # indices of free accelerations

  # Initial guesses
  # x0 = Vector{Float64}(j[1].RB1.x) # Initial guess is provided initial location
  x0 = []; # Initial guess
  for k = 1:length(j)
    append!(x0,j[k].RB2.x)
  end
  extFListCopy = externalFTotal(0.0,j);
  u0 = genU(extFListCopy); # Initial guess for control input

  # n_xu = Integer(length(x0) + length(u0)/2); # Number of x,u (no. of optimization variables excluding gamma)
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
  function _nonlinear_constraints(z::AbstractArray{T}) where T<:Real
    xVar, uVar = getXU(z)
    x = [Vector{T}(j[1].RB1.x); xVar]
    u = [zeros(T,6) reshape(uVar,(6, nB-1))]
    # x = [Vector{T}(j[1].RB1.x);z]; # First body is always the inertial frame
    # x = z[1:14*(nB-1)];
    # u2 = reshape(z[14*(nB-1)+1:end-1],(6,nB-1));
    # u = [zeros(T,6) u2];
    # u = zeros(6,nB); # No external forces or torques in equilibrium
    xdot = fxdot(x,u,j,GravityInInertial)

    ## Relevant quantities (second derivatives only - zero acceleration)
    y = Vector{T}(undef,7*nB);
    for i=1:nB
      y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
    end
    deleteat!(y,1:7); # inertial frame has no acceleration
    deleteat!(y,freeAccInd);

    # Quatnorm constraint
    q_constr = zeros(T,nB-1)
    for i=1:nB-1
      q_constr[i] = norm(xVar[14*(i-1)+4:14*(i-1)+7])^2 - 1;
    end

    # JointLocation constraint  (translation constraint)
    jLoc_constr = zeros(T,3,length(j))
    for i=1:length(j)
      b1_id = j[i].RB1.bodyID; b2_id = j[i].RB2.bodyID;
      x1 = x[(b1_id-1)*14+1:b1_id*14]
      x2 = x[(b2_id-1)*14+1:b2_id*14]
      j_b1 = x1[1:3] + transpose(quat2dcm(x1[4:7]))*j[i].pos1;
      j_b2 = x2[1:3] + transpose(quat2dcm(x2[4:7]))*j[i].pos2;
      jLoc_constr[:,i] = j_b1 - j_b2;
    end

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
    # cfg1 = ForwardDiff.GradientConfig(eval_f,z,ForwardDiff.Chunk{1}());
    grad_f[:] = ForwardDiff.gradient(eval_f,z)#;, cfg1)
    # return grad_f
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
        cfg1 = ForwardDiff.GradientConfig(_nonlinear_constraints, z, ForwardDiff.Chunk{1}());
        values[:] = ForwardDiff.gradient(_nonlinear_constraints,z, cfg1)
    end

    # return jac_g
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
      cfg2 = ForwardDiff.HessianConfig(_nonlinear_constraints, z, ForwardDiff.Chunk{2}());
      H = obj_factor*ForwardDiff.hessian(eval_f,z);
      for i=1:length(lambda)
        H += lambda[i]*ForwardDiff.hessian(_nonlinear_constraints, z, cfg2);
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

  # n = Integer(14*(nB-1) + 1); # number of optimization variables (discounting inertial frame)
  # inGuess = [x0;0.0];
  # n = Integer(20*(nB-1) + 1); # number of optimization variables
  # inGuess = [x0;u0[:,2:end][:];0.0];
  nele_jac = Integer(m*n);
  nele_hess = Integer(n*(n+1)/2);
  prob = createProblem(n, lb, ub, m, lb_g, ub_g, nele_jac, nele_hess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);
  prob.x = inGuess; # Initial guess.
  status = solveProblem(prob);

  trim_x = prob.x[1:length(x0)];
  gam = prob.x[end];
  trim_u = prob.x[length(x0)+1:end-1];
  return trim_x, trim_u, gam
  # return trim_x, gam
end
##
ix= zeros(Integer,20*(nB-1)+1);
## WITH U=0
# for i=1:nB-1
#   ix[20*(i-1)+15:20*i] .= 1
# end
ix[15:20] .= 1;
# ix[29:40] .= 1;
trim_U = zeros(6,nB)
##
iy = zeros(Integer,7);
out = trimDiff(j,g,ix=ix, iy=iy);
# trim_x, gam = out;
trim_x, trim_u, gam = out;
println("trim_x = ", trim_x)
trim_X = [j[1].RB1.x;trim_x]
# trim_U = [zeros(6) reshape(trim_u,(6,nB-1))];
# println("trim_U = ", trim_U[:,2])
# trim_U = zeros(6,nB);
println("qconstr = ", norm(trim_x[4:7])-1)
println("ẍ =", norm(fxdot(trim_X,trim_U, j, g)))
##
# updateRigidBody!(j[1].RB2,trim_x)
##
# function _nonlinear_constraints(z::AbstractArray{T}) where T<:Real
#   nB = Int64((length(z)-1)/20); # Number of bodies
#   x = z[1:14*nB]
#   u = reshape(z[14*nB+1:end-1],(6,nB));
#   xdot = fxdot(x,u,j,GravityInInertial)
#   # Relevant quantities (second derivatives only - zero acceleration)
#   y = Vector{T}(undef,7*nB);
#   y[1:7] = xdot[8:14];
#   for i=1:nB
#     y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
#   end
#
#   J = sum(y.*y) - z[end]
#   return J
# end
#
# temp_x,temp_u = getXU_0(j)
# temp_z = [temp_x;temp_u[:];0.0];
##
# temp_constr_z = _nonlinear_constraints(temp_z);
# println("temp_constr_z = ", temp_constr_z);
# grad_nc = z->ForwardDiff.gradient(_nonlinear_constraints,z);
# @time grad_nc(temp_z);
# hess_nc = z->ForwardDiff.hessian(_nonlinear_constraints,z);
# hess_nc(temp_z)
# cfg1 = ForwardDiff.HessianConfig(_nonlinear_constraints,temp_z,ForwardDiff.Chunk{10}());
# @time ForwardDiff.hessian(_nonlinear_constraints, temp_z, cfg1);
##
# function tempFn(z)
#   # return sqrt(diagm(z))
#   return (diagm(z))^(0.5)
# end
#
# temp_z = rand(5)
# # println(norm(sqrt(temp_z) - temp_z^(0.5)))
# # tempFn(temp_z)
# ForwardDiff.jacobian(tempFn,temp_z)
##
# clearconsole();
# function tempFn(z)
#   return transpose(z)*z
# end
#
# h = z->ForwardDiff.hessian(tempFn,z);
# # h2z) = z->ForwardDiff.hessian(tempFn,z);
# temp_z2 = rand(40);
# @time h(temp_z2)
# @time ForwardDiff.hessian(tempFn,temp_z2)
