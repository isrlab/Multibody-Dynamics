## test DoublePendulum after creating a joint between two single pendulums
# include("../src/plotSol.jl")
include("../src/plotSol.jl")
include("../src/simulateDiff.jl")
include("../src/OrientationConversion.jl")
using Revise
using JLD
using BenchmarkTools
##
clearconsole()
j1 = Nothing
j2 = Nothing
## Pendulum 1
m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1]

# Testing Dynamics with Revolute Joint
Pend1 = RigidBody(m,I1,2)
RbI = InertialFrameAsRB()
x0R1 = [l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Pend1,x0R1)

axisY = [0.0 1.0 0.0][:] # Axis about which bar is revolving

rj1 = [0.0 0.0 0.0][:] # Joint Location in body frame of first body
rj2 = [-l/2 0.0 0.0][:] # Joint location in body frame of second body

j1_1 = Joint(RbI,Pend1,rj1,rj2,type="Revolute",axis=axisY)
j1 = [j1_1];

# just simulate without any external force or gravity first
g = zeros(3);
tEnd = 1.0; tInt = 0.01;
tSim, solFinal = simulateDiff(tEnd,tInt,j1,g=g);
solPend1 = solFinal[1];

## Pendulum 2
Pend2 = RigidBody(m,I1,2);
x0R2 = [3*l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Pend2,x0R2)

rj1 = [l 0.0 0.0][:] #Joint location wrt inertial frame
rj2 = [-l/2 0.0 0.0][:] #Joint location in body frame of pend2

j2_1 = Joint(RbI,Pend2,rj1,rj2,type="Revolute",axis=axisY);
j2 = [j2_1];

# just simulate without any external force or gravity first
g = zeros(3);
tEnd = 1.0; tInt = 0.01;
tSim, solFinal = simulateDiff(tEnd,tInt,j1,g = g);
solPend2 = solFinal[1];

## Join the two pendulums
println("########################")
println("Joining the two pendulums")
rj1 = [l/2;zeros(2)];
rj2 = [-l/2;zeros(2)];
joinBodies!(j1,j2,2,2,rj1,rj2, type = "Revolute", axis = axisY);

## Simulation
g = [0.0,0.0,-9.806]; # Gravity Vector.
tEnd = 1.0
tSpan = 0.01
tSim, solFinal = simulateDiff(tEnd,tSpan,j1,g=g)

solPend1  = solFinal[1]
solPend2  = solFinal[2]
##
# Check if joint location has moved
jointLoc1 = Matrix{Float64}(undef,length(tSim),3)
jointLoc2 = Matrix{Float64}(undef,length(tSim),3)
ωSol = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    Pend1.dcm = quat2dcm(solPend1.β[i,:])
    Pend2.dcm = quat2dcm(solPend2.β[i,:])
    # ωSol[i,:] = angVel(solPend1.β[i,:],solPend1.βdot[i,:])
    jointLoc1[i,:] = solPend1.r[i,:] + transpose(Pend1.dcm)*rj1
    jointLoc2[i,:] = solPend2.r[i,:] + transpose(Pend2.dcm)*rj2
end
##
close("all");
# plotErrNorm(tSim,solPend.β)
plotPos(tSim,jointLoc1 .- jointLoc2)
# plotPos(tSim,jointLoc2)
# plotQuat(tSim, solPend1.β)
# plotAngVel(tSim, ωSol)
# plotPos(tSim, solPend.r)
## Compare with MinReal
solMin = load("DoublePendulumMin.jld")
r1Min = solMin["r1Min"]
r2Min = solMin["r2Min"]

# Plotting errors
plotPos(tSim,r1Min-solPend1.r)
plotPos(tSim,r2Min-solPend2.r)
# plotAngVel(tSim,ωMin+ωSol)
# plotPos(tSim,solPend.r)
