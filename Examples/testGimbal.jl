# To test a simplified 2 axis gimbal.
# Modeling the system as a cube and the gimbal placed on top.
# We need two revolute joints between the cube and gimbal.

include("../src/plotSol.jl")
include("../src/simulate.jl")
include("../src/OrientationConversion.jl")
clearconsole()

## Inertial Properties
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

## Initialisation
InFrame = InertialFrameAsRB()
Body = RigidBody(mb,Ib,2)
Gimbal = RigidBody(mg,Ig,3)

x0Body = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Body,x0Body)
x0Gimbal = [[0.0;0.0;0.143];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Gimbal,x0Gimbal)

## Joint Descriptions
# Joint 1 (Free) between the inertial frame and body.
j1 = Joint(InFrame,Body,zeros(3),zeros(3))

# Joint 2 (Revolute2) between the body and gimbal
axisZ = [0.0 0.0 1.0][:]
rjBody = [0.0 0.0 0.134][:]
rjGimbal = [0.0 0.0 -0.009][:]
j2 = Joint(Body,Gimbal,rjBody,rjGimbal,
     type = "Revolute",axis = axisZ, jointTorque = [0.0 0.0 0.0][:])

## Simulation
tEnd = 1.0
tSpan = 0.01
g = [0.0;0.0;0.0]
tSim, solFinal = simulate(tEnd,tSpan,j1,j2,g=g)#,extFVec = extFList)

solBody = solFinal[1]
solGimbal = solFinal[2]

## Plotting
plotErrNorm(tSim,solBody.β)
plotErrNorm(tSim,solGimbal.β)
# Check if joint location has moved
jointLoc = Matrix{Float64}(undef,length(tSim),3)
ωBody = Matrix{Float64}(undef,length(tSim),3)
ωGimbal = Matrix{Float64}(undef,length(tSim),3)
ωGimbalInBody = Matrix{Float64}(undef,length(tSim),3)
ωRel = Matrix{Float64}(undef,length(tSim),3) # Relative
for i=1:length(tSim)
    Body.dcm = quat2dcm(solBody.β[i,:])
    Gimbal.dcm = quat2dcm(solGimbal.β[i,:])
    ωBody[i,:] = angVel(solBody.β[i,:],solBody.βdot[i,:])
    ωGimbal[i,:] = angVel(solGimbal.β[i,:],solGimbal.βdot[i,:])
    ωGimbalInBody[i,:] = quat2dcm(solBody.β[i,:])*
                         transpose(quat2dcm(solGimbal.β[i,:]))*
                         ωGimbal[i,:]
    ωRel[i,:] = ωGimbalInBody[i,:] - ωBody[i,:]
    jointLoc[i,:] = solBody.r[i,:] + transpose(Body.dcm)*rjBody -
                    solGimbal.r[i,:] - transpose(Gimbal.dcm)*rjGimbal
end
plotPos(tSim,jointLoc)
plotPos(tSim,solBody.r)
plotPos(tSim,solGimbal.r)
plotVel(tSim,solBody.v)
plotVel(tSim,solGimbal.v)
plotAngVel(tSim,ωBody)
plotAngVel(tSim,ωGimbal)
plotAngVel(tSim,ωGimbalInBody)
plotAngVel(tSim,ωRel)
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
