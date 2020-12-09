# include("../src/plotSol.jl")
include("../src/plotSol.jl")
include("../src/simulate.jl")
include("../src/OrientationConversion.jl")
using Revise
using JLD
using BenchmarkTools
clearconsole()
##
m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1]

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,2)
R2 = RigidBody(m,I1,3)
RbI = InertialFrameAsRB()
x0R1 = [l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(R1,x0R1)
x0R2 = [3*l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(R2,x0R2)

axisY = [0.0 1.0 0.0][:] # Axis about which bar is revolving

rj1 = [0.0 0.0 0.0][:] # Joint Location in body frame of first body
rj2 = [-l/2 0.0 0.0][:] # Joint location in body frame of second body

rj3 = [l/2 0.0 0.0][:] #Joint location in body frame of pend1
rj4 = [-l/2 0.0 0.0][:] #Joint location in body frame of pend2

j1 = Joint(RbI,R1,rj1,rj2,type="Revolute",axis=axisY)
j2 = Joint(R1,R2,rj3,rj4,type="Revolute",axis=axisY)

# External Forces Definition
g = MVector{3}([0.0,0.0,-9.806]) # Gravity Vector.

## Simulation
tEnd = 1.0
tSpan = 0.01
tSim, solFinal = simulate(tEnd,tSpan,j1,j2,g=g)

solPend1  = solFinal[1]
solPend2  = solFinal[2]
##
# ωCube, ωProp1, ωRel1, jointLoc1 = checkRevJoint(solQuad, solProp1, rjCube1, rjProp);
# ωCube, ωProp2, ωRel2, jointLoc2 = checkRevJoint(solQuad, solProp2, rjCube2, rjProp);
##
# After adding solSim

# Check if joint location has moved
jointLoc1 = Matrix{Float64}(undef,length(tSim),3)
jointLoc2 = Matrix{Float64}(undef,length(tSim),3)
ωSol = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    R1.dcm = quat2dcm(solPend1.β[i,:])
    R2.dcm = quat2dcm(solPend2.β[i,:])
    # ωSol[i,:] = angVel(solPend1.β[i,:],solPend1.βdot[i,:])
    jointLoc1[i,:] = solPend1.r[i,:] + transpose(R1.dcm)*rj3
    jointLoc2[i,:] = solPend2.r[i,:] + transpose(R2.dcm)*rj4
end
##
close("all");
# plotErrNorm(tSim,solPend.β)
plotPos(tSim,jointLoc1)# .- jointLoc2)
plotPos(tSim,jointLoc2)
# plotQuat(tSim, solPend1.β)
# plotAngVel(tSim, ωSol)
# plotPos(tSim, solPend.r)
## Compare with MinReal
solMin = load("PendulumMin.jld")
ωMin = solMin["ωMin"]
rMin = solMin["rMin"]
vMin = solMin["vMin"]

# Plotting errors
plotPos(tSim,rMin-solFinal[1].r)
plotVel(tSim,vMin-solFinal[1].v)
plotAngVel(tSim,ωMin+ωSol)
plotPos(tSim,solPend.r)

##
j1 = Nothing
j2 = Nothing

#
# # tSim = sol.t
# # rSol = transpose(sol[15:17,:])
# # vSol = transpose(sol[22:24,:])
# # βSol = transpose(sol[18:21,:])
# # βdotSol = transpose(sol[25:28,:]);
#
# # Check errNorm
# plotErrNorm(tSim,βSol)
# # Check if joint location has moved
# jointLoc = Matrix{Float64}(undef,length(tSim),3)
# ωSol = Matrix{Float64}(undef,length(tSim),3)
# for i=1:length(tSim)
#     R1.dcm = quat2dcm(βSol[i,:])
#     ωSol[i,:] = angVel(βSol[i,:],βdotSol[i,:])
#     jointLoc[i,:] = rSol[i,:] + transpose(R1.dcm)*rj2
# end
# plotPos(tSim,jointLoc)
#
# # Compare with MinReal
# solMin = load("PendulumMin.jld")
# ωMin = solMin["ωMin"]
# rMin = solMin["rMin"]
# vMin = solMin["vMin"]
#
# # Plotting errors
# plotPos(tSim,rMin-rSol)
# plotVel(tSim,vMin-vSol)
# plotAngVel(tSim,ωMin+ωSol)
