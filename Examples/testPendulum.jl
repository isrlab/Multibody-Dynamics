include("../src/plotSol.jl")
include("../src/simulateDiff.jl")
include("../src/OrientationConversion.jl")
using Revise
using JLD
using BenchmarkTools
# using Debugger
##
clearconsole()
j1 = Nothing
##
m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1];

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,2)
RbI = InertialFrameAsRB()

# Suspended horizontally
x0R1 = ([l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]);
initialiseRigidBody!(R1,x0R1)

# Equilibrium point
# x0R1 = ([zeros(2);-l/2;[1;zeros(3)];zeros(3);zeros(4)]);
# initialiseRigidBody!(R1,x0R1)

# Axis about which bar is revolving
axisY = [0.0 1.0 0.0][:];


rj1 = [0.0 0.0 0.0][:]; # Joint Location in body frame of first body
rj2 = [-l/2 0.0 0.0][:]; # Joint location in body frame of second body

j1 = Joint(RbI,R1,rj1,rj2,type="Revolute",axis=axisY)

# External Forces Definition
g = [0.0,0.0,-9.806]; # Gravity Vector.

j = [j1];


## Linearization
x0, u0 = getXU_0(j);
A,B = linearizeDiff(x0, u0, j);
println("A = ",A)
println("B = ",B)
# sleep(1000);

## Simulation
tEnd = (1.0);
tInt = (0.01); # Only for fixed step solver, currently using adaptive step.

@btime simulateDiff(tEnd, tInt, j, g=g);
sleep(1000);
tSim, solFinal = simulateDiff(tEnd,tInt,j,g=g)

solPend1  = solFinal[1];
##
# Check if joint location has moved
ωSol, jointLoc1 = checkRevJointIn(solPend1, rj2)
# jointLoc1 = Matrix{Float64}(undef,length(tSim),3)
# ωSol = Matrix{Float64}(undef,length(tSim),3)
# for i=1:length(tSim)

#     R1.dcm = quat2dcm(solPend1.β[i,:])
#     ωSol[i,:] = angVel(solPend1.β[i,:],solPend1.βdot[i,:])
#     jointLoc1[i,:] = solPend1.r[i,:] + transpose(R1.dcm)*rj2
# end
## Plotting Errors
close("all");

plotErrNorm(tSim,solPend1.β) # Drift in uaternion unit norm constraint
plotPos(tSim,jointLoc1) # Location of joint connecting inertial frame to pendulum. Should be zero.

# plotQuat(tSim, solPend1.β)
plotAngVel(tSim, ωSol)
plotPos(tSim, solPend1.r)
## Compare with MinReal
# Can only be done if simulated with fixed-step
# solMin = load("PendulumMin.jld")
# ωMin = solMin["ωMin"]
# rMin = solMin["rMin"]
# vMin = solMin["vMin"]
#
# # Plotting errors
# plotPos(tSim,rMin-solFinal[1].r)
# plotVel(tSim,vMin-solFinal[1].v)
# plotAngVel(tSim,ωMin+ωSol)
# plotPos(tSim,solPend1.r)
