include("plotSol.jl")
include("simulate.jl")
include("OrientationConversion.jl")
using Revise
using JLD
clearconsole()

m = 1.0; l = 1.0; r = 0.01 # Mass, length, and radius of bar
# Assuming bar revolves about Z axis
I1 = [m*r^2/2 0 0; 0 m*(l^2/12 + r^2/4) 0; 0 0 m*(l^2/12 + r^2/4)]

# Testing Dynamics with Revolute Joint
# System looks like this:
# ====|====     | is the joint axis.
RbI = InertialFrameAsRB()
R1 = RigidBody(m,I1,"quaternions")
x0R1 = [l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
R1 = initialiseRigidBody(R1,x0R1)

axis = [0.0 0.0 1.0][:] # Axis about which bar is revolving

rj1 = [0.0 0.0 0.0][:] # Joint Location in body frame of first body
rj2 = [0 0.0 0.0][:] # Joint location in body frame of second body

j = Joint(RbI,R1,rj1,rj2,type="Revolute",axis=axis)

# Simulation
tEnd = 5.0
tSpan = 0.01
g = [0.0;0.0;-9.806]
tSim, solFinal = simulate(tEnd,tSpan,j,g=g)

solPend  = solFinal[1]

# After adding solSim
plotErrNorm(tSim,solPend.β)
# Check if joint location has moved
jointLoc = Matrix{Float64}(undef,length(tSim),3)
ωSol = Matrix{Float64}(undef,length(tSim),3)
rEnd = Matrix{Float64}(undef,length(tSim),3)
endPos = [l/2;0;0]
for i=1:length(tSim)
    R1.dcm = quat2dcm(solFinal[1].β[i,:])
    ωSol[i,:] = angVel(solFinal[1].β[i,:],solFinal[1].βdot[i,:])
    rEnd[i,:] = solPend.r[i,:] + transpose(R1.dcm)*endPos
    jointLoc[i,:] = solFinal[1].r[i,:] + transpose(R1.dcm)*rj2
end
plotPos(tSim,jointLoc)
plotPos(tSim,solPend.r)
plotAngVel(tSim,ωSol)
plotPos(tSim,rEnd)
plotQuat(tSim,solPend.β)
# plotQuatDot(tSim,solPend.βdot)
