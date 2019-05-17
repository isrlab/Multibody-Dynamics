# To test a simplified co ax copter.
# Modeling the system as a cube, gimbal and 2 propellers.
# A revolute2 joint between the cube and gimbal
# 1 Rev Joint each between the gimbal and either propeller.

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
Prop1 = RigidBody(mp,Ip,4)
Prop2 = RigidBody(mp,Ip,5)

x0Body = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
Body = initialiseRigidBody(Body,x0Body)

x0Gimbal = [[0.0;0.0;0.143];[1;zeros(3)];zeros(3);zeros(4)]
Gimbal = initialiseRigidBody(Gimbal,x0Gimbal)

x0Prop1 = [[0.0;0.0;0.234];[1;zeros(3)];zeros(3);zeros(4)]
Prop1 = initialiseRigidBody(Prop1,x0Prop1)

x0Prop2 = [[0.0;0.0;0.24];[1;zeros(3)];zeros(3);zeros(4)]
Prop2 = initialiseRigidBody(Prop2,x0Prop2)

## Joint Descriptions
# Joint 1 (Free) between the inertial frame and body.
j1 = Joint(InFrame,Body,zeros(3),zeros(3))

# Joint 2 (Revolute2) between the body and gimbal
axisZ = [0.0 0.0 1.0][:]
rjBody = [0.0 0.0 0.134][:]
rjGimbal1 = [0.0 0.0 -0.009][:]
j2 = Joint(Body,Gimbal,rjBody,rjGimbal1,
     type = "Revolute",axis = axisZ)#, jointTorque = [0.00 0.0 0.00][:])

# Joint 2 (Revolute2) between the body and gimbal
axisZ = [0.0 0.0 1.0][:]
rjGimbal2 = [0.0 0.0 0.091][:]
rjProp1 = [0.0 0.0 0.0][:]
j3 = Joint(Gimbal,Prop1,rjGimbal2,rjProp1,
     type = "Revolute",axis = axisZ, jointTorque = [0.0 0.0 0.00][:])


# Joint 2 (Revolute2) between the body and gimbal
axisZ = [0.0 0.0 1.0][:]
rjGimbal3 = [0.0 0.0 0.097][:]
rjProp2 = [0.0 0.0 0.0][:]
j4 = Joint(Gimbal,Prop2,rjGimbal3,rjProp2,
     type = "Revolute", axis = axisZ, jointTorque = [0.0 0.0 0.00][:])


## Simulation
tEnd = 0.1
tSpan = 0.01
g = [0.0;0.0;0.0]
tSim, solFinal = simulate(tEnd,tSpan,j1,j2,j3,j4,g=g)#,extFVec = extFList)

solBody = solFinal[1]
solGimbal = solFinal[2]
solProp1 = solFinal[3]
solProp2 = solFinal[4]

## Check constraints
# Check if joint location has moved
jointLoc2 = Matrix{Float64}(undef,length(tSim),3)
jointLoc3 = Matrix{Float64}(undef,length(tSim),3)
jointLoc4 = Matrix{Float64}(undef,length(tSim),3)
ωBody = Matrix{Float64}(undef,length(tSim),3)
ωGimbal = Matrix{Float64}(undef,length(tSim),3)
ωProp1 = Matrix{Float64}(undef,length(tSim),3)
ωProp2 = Matrix{Float64}(undef,length(tSim),3)
ωGimbalInBody = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    Body.dcm = quat2dcm(solBody.β[i,:])
    Gimbal.dcm = quat2dcm(solGimbal.β[i,:])
    Prop1.dcm = quat2dcm(solProp1.β[i,:])
    Prop2.dcm = quat2dcm(solProp2.β[i,:])
    ωBody[i,:] = angVel(solBody.β[i,:],solBody.βdot[i,:])
    ωGimbal[i,:] = angVel(solGimbal.β[i,:],solGimbal.βdot[i,:])
    ωProp1[i,:] = angVel(solProp1.β[i,:],solProp1.βdot[i,:])
    ωProp2[i,:] = angVel(solProp2.β[i,:],solProp2.βdot[i,:])
    ωGimbalInBody[i,:] = quat2dcm(solBody.β[i,:])*
                         transpose(quat2dcm(solGimbal.β[i,:]))*
                         ωGimbal[i,:]
    jointLoc2[i,:] = solBody.r[i,:] + transpose(Body.dcm)*rjBody -
                     solGimbal.r[i,:] -
                     transpose(Gimbal.dcm)*rjGimbal1
    jointLoc3[i,:] = solGimbal.r[i,:] +
                     transpose(Gimbal.dcm)*rjGimbal2 -
                     solProp1.r[i,:] - transpose(Prop1.dcm)*rjProp1
    jointLoc4[i,:] = solGimbal.r[i,:] +
                     transpose(Gimbal.dcm)*rjGimbal3 -
                     solProp2.r[i,:] - transpose(Prop2.dcm)*rjProp2
end

## Plotting
# Joint Locations
plotPos(tSim,jointLoc2)
plotPos(tSim,solBody.r)
plotPos(tSim,solGimbal.r)

plotErrNorm(tSim,solBody.β)
plotErrNorm(tSim,solGimbal.β)
plotErrNorm(tSim,solProp1.β)
plotErrNorm(tSim,solProp2.β)

plotAngVel(tSim,ωBody)
plotAngVel(tSim,ωGimbal)
# plotAngVel(tSim,ωProp1)
# plotAngVel(tSim,ωProp2)
# plotAngVel(tSim,ωGimbalInBody)
# plotAngVel(tSim,ωCube - ωProp)
