# To test a simplified quadrotor.
# Modeling the quadrotor frame as a cube and the props
# as another rigid body with a revolute joint between the two.

include("../src/plotSol.jl")
include("../src/simulate.jl")
include("../src/OrientationConversion.jl")
clearconsole()
mb = 0.155; mp = 0.02;

Ib =  [ 0.0009    0.0000   -0.0000;
        0.0000    0.0053    0.0000;
       -0.0000    0.0000    0.0054]#*(10^-3)

Ip =    [ 0.0029    0.0000    0.0000;
          0.0000    0.5824    0.0000;
          0.0000    0.0000    0.5805]*(10^-4)

InFrame = InertialFrameAsRB()
QuadCube = RigidBody(mb,Ib,"quaternions")
Props = RigidBody(mp,Ip,"quaternions")

x0Cube = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
QuadCube = initialiseRigidBody(QuadCube,x0Cube)
x0Props = [[0.0;0.0;0.143];[1;zeros(3)];zeros(3);zeros(4)]
Props = initialiseRigidBody(Props,x0Props)

axis = [0.0 0.0 1.0][:]
rjCube = [0.0 0.0 0.143][:]
rjProp = [0.0 0.0 0.0][:]

j1 = Joint(InFrame,QuadCube,zeros(3),zeros(3))
j2 = Joint(QuadCube,Props,rjCube,rjProp,type = "Revolute",axis = axis,jointTorque = [0.0;0.0;0.0001])

# External Forces Definition
g = [0.0;0.0;0.0]
extFList = Vector{extForces}(undef,3)
extFList[1] = zeroExtForce()
extFList[2] = extForces(transpose((j1.RB2.m)*-g),zeros(1,3),[0.0 0.0 0.0])
extFList[3] = extForces(transpose((j2.RB2.m)*-g),zeros(1,3),[0.0 0.0 0.0])

# Simulation
tEnd = 0.1
tSpan = 0.01
g = [0.0;0.0;0.0]
tSim, solFinal = simulate(tEnd,tSpan,j1,j2,g=g,extFVec = extFList)

solQuad = solFinal[1]
solProp = solFinal[2]

plotErrNorm(tSim,solQuad.β)
plotErrNorm(tSim,solProp.β)
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
