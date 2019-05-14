# To test a simplified 2 axis gimbal.
# Modeling the system as a cube and the gimbal placed on top.
# We need two revolute joints between the cube and gimbal.

include("../src/plotSol.jl")
include("../src/simulate.jl")
include("../src/OrientationConversion.jl")
clearconsole()
mb = 0.755; mg = 0.215; mp = 0.02;

Ib =  [ 0.0009    0.0000   -0.0000;
        0.0000    0.0053    0.0000;
       -0.0000    0.0000    0.0054]#*(10^-3)

Ip =    [ 0.0029    0.0000    0.0000;
          0.0000    0.5824    0.0000;
          0.0000    0.0000    0.5805]*(10^-4)

Ig =    [ 0.1435    0.0000    0.0000;
          0.0000    0.0498    0.0000;
          0.0000    0.0000    0.1381]*(10^-3)

InFrame = InertialFrameAsRB()
Body = RigidBody(mb,Ib,"quaternions")
Gimbal = RigidBody(mg,Ig,"quaternions")
VirtualLink = RigidBody(0.0,zeros(3,3),"quaternions")

x0Body = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
Body = initialiseRigidBody(Body,x0Body)
x0Gimbal = [[0.0;0.0;0.143];[1;zeros(3)];zeros(3);zeros(4)]
Gimbal = initialiseRigidBody(Gimbal,x0Gimbal)
VirtualLink = initialiseRigidBody(VirtualLink,x0Gimbal)

# Joint 1 (Free) between the inertial frame and body.
j1 = Joint(InFrame,QuadCube,zeros(3),zeros(3))

# Joint 2
axisX = [1.0 0.0 0.0][:]
rjBody = [0.0 0.0 0.143][:]
rjVirtualLink = [0.0 0.0 0.0][:]
j2 = Joint(Body,VirtualLink,rjBody,rjVirtualLink,
     type = "Revolute",axis = axisX,jointTorque = [0.0 0.0 0.0][:])

# Joint 3
axisY = [0.0 1.0 0.0][:]
rjVirtualLink = [0.0 0.0 0.0][:]
rjGimbal = [0.0 0.0 0.0][:]
j3 = Joint(VirtualLink,Gimbal,rjVirtualLink,rjGimbal,
     type="Revolute",axis = axisY,jointTorque = [0.0 0.0 0.0][:])


# External Forces Definition
g = [0.0;0.0;0.0]
extFList = Vector{extForces}(undef,4)
extFList[1] = zeroExtForce()
extFList[2] = zeroExtForce()
extFList[3] = zeroExtForce()
extFList[4] = zeroExtForce()

# Simulation
tEnd = 1.0
tSpan = 0.01
g = [0.0;0.0;0.0]
tSim, solFinal = simulate(tEnd,tSpan,j1,j2,g=g,extFVec = extFList)

solBody = solFinal[1]
solVirtualLink = solFinal[2]
solGimbal = solFinal[3]

plotErrNorm(tSim,solBody.β)
plotErrNorm(tSim,solGimbal.β)
# Check if joint location has moved
jointLoc = Matrix{Float64}(undef,length(tSim),3)
ωCube = Matrix{Float64}(undef,length(tSim),3)
ωProp = Matrix{Float64}(undef,length(tSim),3)
ωPropInCube = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    QuadCube.dcm = quat2dcm(solQuad.β[i,:])
    ωCube[i,:] = angVel(solQuad.β[i,:],solQuad.βdot[i,:])
    ωProp[i,:] = angVel(solProp.β[i,:],solProp.βdot[i,:])
    ωPropInCube[i,:] = quat2dcm(solQuad.β[i,:])*transpose(quat2dcm(solProp.β[i,:]))*ωProp[i,:]
    jointLoc[i,:] = solQuad.r[i,:] + transpose(QuadCube.dcm)*rjCube - solProp.r[i,:]
end
plotPos(tSim,jointLoc)
plotVel(tSim,solQuad.v - solProp.v)
plotAngVel(tSim,ωCube)
plotAngVel(tSim,ωProp)
plotAngVel(tSim,ωPropInCube)
# plotAngVel(tSim,ωCube - ωProp)

# # Check revJoint axis
# axisSol = Matrix{Float64}(undef,length(tSim),3)
# errAxis = Matrix{Float64}(undef,length(tSim),3)
# errAxisNorm = Vector{Float64}(undef,length(tSim))
# for i=1:length(tSim)
#     axisSol[i,:] = solProp.β[i,2:4]/norm(solProp.β[i,2:4],2)
#
#     # β1inv = [solQuad.β[i,1];solQuad.β[i,2:4]]
#     # βRev = quaternionProduct(solProp.β[i,:],β1inv)
#     # axisSol[i,:] = βRev[2:4]/norm(βRev[2:4])
#
#     # axisSol[i,:] = (solQuad.β[i,:])*solProp.β[i,2:4]/norm(solProp.β[i,2:4],2)
#
#     # axisSol[i,:] = axixRev(solQuad.β[i,:],solProp.β[i,:])
#     errAxis[i,:] = axisSol[i,:] - axis
#     errAxisNorm[i] = norm(axisSol[i,:],2) - 1.0
# end
# plot(tSim[2:end],errAxisNorm[2:end])
# plotPos(tSim,errAxis)
