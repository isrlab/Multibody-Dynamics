# include("../src/plotSol.jl")
include("../src/simulateDiff.jl")
include("../src/OrientationConversion.jl")
using Revise
# using JLD
using BenchmarkTools
# using Debugger

# clearconsole()
j1 = Nothing

m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1];

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,2)
RbI = InertialFrameAsRB()

# Suspended horizontally
x0R1 = ([l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]);
initialiseRigidBody!(R1,x0R1)

# Perturbed from equilibrium position
# x0R1 = ([[0.3;0.0;-0.4];[1;zeros(3)];zeros(3);zeros(4)]);
# initialiseRigidBody!(R1,x0R1)

# Equilibrium point
# x0R1 = ([zeros(2);-l/2;[1;zeros(3)];zeros(3);zeros(4)]);
# initialiseRigidBody!(R1,x0R1)

# Axis about which bar is revolving
axisY = [0.0 1.0 0.0][:];

rj1 = [0.0 0.0 0.0][:]; # Joint Location in body frame of first body
rj2 = -R1.x[1:3];

j1 = Joint(RbI,R1,rj1,rj2,type="Revolute",axis=axisY)

# External Forces Definition
g = [0.0,0.0,-9.806]; # Gravity Vector.

j = [j1];
##
using Ipopt
using BenchmarkTools
using GenericLinearAlgebra
# clearconsole();
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
function trimDiff(j::Vector{Joint},GravityInInertial::Vector{Float64})
  # Initial guesses
  # x0 = Vector{Float64}(j[1].RB1.x) # Initial guess is provided initial location
  x0 = []; # Initial guess
  for k = 1:length(j)
    append!(x0,j[k].RB2.x)
  end
  extFListCopy = externalFTotal(0.0,j);
  u0 = genU(extFListCopy); # Initial guess for control input

  n_xu = Integer(length(x0) + length(u0)/2); # Number of x,u (no. of optimization variables excluding gamma)
  lb = [-Inf*ones(n_xu);0.0]; ub = [Inf*ones(n_xu);Inf] # Bounds on optimization variables
  # lb[4] = 0.0;
  # lb[14+4] = 0.0;

  m = 1; # Number of nonlinear constraints (xdot = f(x,u)) for second order variables
  lb_g = [-Inf]; ub_g = [0.0];

  function _nonlinear_constraints(z::AbstractArray{T}) where T<:Real
    x = [Vector{T}(j[1].RB1.x);z]; # First body is always the inertial frame
    # x = z[1:14*(nB-1)];
    u2 = reshape(z[14*(nB-1)+1:end-1],(6,nB-1));
    u = [zeros(T,6) u2];
    # u = zeros(6,nB); # No external forces or torques in equilibrium
    xdot = fxdot(x,u,j,GravityInInertial)

    ## Relevant quantities (second derivatives only - zero acceleration)
    y = Vector{T}(undef,7*nB);
    for i=1:nB
      y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
    end
    deleteat!(y,1:7);

    # Quatnorm constraint (currently specific to system with 2 RBs only)
    q_constr1 = norm(z[4:7])^2 - 1;
    # q_constr2 = norm(z[14+4:14+7])^2 - 1;
    q_constr = [q_constr1];#; q_constr2];

    # JointLocation constraint  (currently specific to system with 2 RBs only)
    # j_b1 = z[1:3] + transpose(quat2dcm(z[4:7]))*j[1].pos1;
    # j_b2 = z[14+1:14+3] + transpose(quat2dcm(z[14+4:14+7]))*j[1].pos2;
    # jLoc_constr = j_b1 - j_b2;

    constr = [y;q_constr;];#jLoc_constr];

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
        # cfg1 = ForwardDiff.GradientConfig(_nonlinear_constraints, z, ForwardDiff.Chunk{1}());
        values[:] = ForwardDiff.gradient(_nonlinear_constraints,z)#, cfg1)
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

  # n = Integer(14*(nB-1) + 1); # number of optimization variables (discounting inertial frame)
  # inGuess = [x0;0.0];
  n = Integer(20*(nB-1) + 1); # number of optimization variables
  inGuess = [x0;u0[:,2][:];0.0];
  nele_jac = Integer(m*n);
  nele_hess = Integer(n*(n+1)/2);
  prob = createProblem(n, lb, ub, m, lb_g, ub_g, nele_jac, nele_hess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);
  # addOption(prob,"linear_solver", "ma77");
  prob.x = inGuess; # Initial guess.
  status = solveProblem(prob);

  trim_x = prob.x[1:length(x0)];
  gam = prob.x[end];
  trim_u = prob.x[length(x0)+1:end-1];
  return trim_x, trim_u, gam
  # return trim_x, gam
end
##
out = trimDiff(j,g);
# trim_x, gam = out;
trim_x, trim_u, gam = out;
println("trim_x = ", trim_x)
trim_X = [j[1].RB1.x;trim_x]
trim_U = [zeros(6) reshape(trim_u,(6,nB-1))];
# println("trim_U = ", trim_U[:,2])
# trim_U = zeros(6,nB);
norm(fxdot(trim_X,trim_U, j, g))
##
# function _nonlinear_constraints(z::Vector{T}) where T<:Real
#   nB = Int64((length(z)-1)/20); # Number of bodies
#   x = z[1:14*nB]
#   u = reshape(z[14*nB+1:end-1],(6,nB));
#   xdot = fxdot(x,u,j,GravityInInertial)
#   # Relevant quantities (second derivatives only)
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
# # clearconsole();
# println("#############")
# println("Computing Hessian")
# cfg1 = ForwardDiff.HessianConfig(_nonlinear_constraints,temp_z,ForwardDiff.Chunk{1}());
# println("With chunk size of 1");
# @btime ForwardDiff.hessian(_nonlinear_constraints, temp_z, cfg1);
# cfg2 = ForwardDiff.HessianConfig(_nonlinear_constraints,temp_z,ForwardDiff.Chunk{2}());
# println("With chunk size of 2");
# @btime ForwardDiff.hessian(_nonlinear_constraints, temp_z, cfg2);
#
# cfg5 = ForwardDiff.HessianConfig(_nonlinear_constraints,temp_z,ForwardDiff.Chunk{5}());
# println("With chunk size of 5");
# @btime ForwardDiff.hessian(_nonlinear_constraints, temp_z, cfg5);
##
# temp_gradz = ForwardDiff.gradient(x -> _nonlinear_constraints(x), temp_z);
# println("temp_gradz = ", temp_gradz)
# temp_hess(z) = ForwardDiff.hessian(_nonlinear_constraints,z)
# @time temp_hess(temp_z)
# @time ForwardDiff.hessian
# temp_hessz = ForwardDiff.hessian(x -> _nonlinear_constraints(x), temp_z);
# println("temp_hessz = ", temp_hessz)
