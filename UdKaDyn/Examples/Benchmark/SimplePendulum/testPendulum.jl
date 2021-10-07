# include("../../../src/plotSol.jl")
# include("../../../src/simulate.jl")
# include("../../../src/linearize.jl")
# include("../../../src/OrientationConversion.jl")
# include("../../../src/trim_kronLazy.jl");
using UdKaDyn, LinearAlgebra
# using Revise, UdKaDyn
# using BenchmarkTools

include("simplePendulumMinReal.jl")
##
# clearconsole()
j1 = Nothing
##
m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m * l^2 / 12 0; 0 0 1];

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m, I1, 2)
RbI = InertialFrameAsRB()

# Suspended at an angle theta from vertical
theta = deg2rad(30)
x_temp = [l/2*sin(theta); 0.0; -l/2*cos(theta)];
x0R1 = ([x_temp;[1;zeros(3)];zeros(3);zeros(4)]);
initialiseRigidBody!(R1,x0R1)

# Axis about which bar is revolving
axisY = [0.0 1.0 0.0][:];

rj1 = [0.0 0.0 0.0][:]; # Joint Location in body frame of first body
rj2 = -R1.x[1:3];

j1 = Joint(RbI, R1, rj1, rj2, type="Revolute", axis=axisY)

# External Forces Definition
g = [0.0,0.0,-9.806]; # Gravity Vector.

j = [j1];

## Linearization
# x0, u0 = getXU_0(j);
# A,B = linearize(x0, u0, j, g);

## Simulation and checking
tEnd = (10.0);
tInt = 1e-2; # Only for fixed step solver, currently using adaptive step.
tSim, solFinal = simulate(tEnd, tInt, j, g=g)
solPend1  = solFinal[1];
println("Simulation complete., plotting results:");

##
# Check if joint location has moved
ωSol, jointLoc1 = checkRevJointIn(solPend1, rj2)
## Plots
close("all");
# errors
plotErrNorm(tSim,solPend1.β) # Drift in quaternion unit norm constraint
plotPos(tSim,jointLoc1) # Location of joint connecting inertial frame to pendulum. Should be zero.

# trajectory
plotQuat(tSim, solPend1.β)
plotPos(tSim, solPend1.r)
plotAngVel(tSim, ωSol)

## Compare with MinReal
rMin, vMin, ωMin = PendulumMinReal(m, l, theta, tInt, tEnd);
# # Plotting errors
plotPos(tSim,rMin-solPend1.r)
plotVel(tSim,vMin-solPend1.v)
plotAngVel(tSim,ωMin-ωSol)
