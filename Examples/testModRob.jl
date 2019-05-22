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
Body1 = RigidBody(mb,Ib,2)
Gimbal1 = RigidBody(mg,Ig,3)
Prop11 = RigidBody(mp,Ip,4)
Prop12 = RigidBody(mp,Ip,5)

Box = RigidBody(mb,Ib,6)

Body2 = RigidBody(mb,Ib,7)
Gimbal2 = RigidBody(mg,Ig,8)
Prop21 = RigidBody(mp,Ip,9)
Prop22 = RigidBody(mp,Ip,10)

# CoAxCop1
x0Body1 = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Body1,x0Body1)

x0Gimbal1 = [[0.0;0.0;0.143];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Gimbal1,x0Gimbal1)

x0Prop11 = [[0.0;0.0;0.234];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop11,x0Prop11)

x0Prop12 = [[0.0;0.0;0.24];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop12,x0Prop12)

# Box
x0Box = [[0.08;0.0;0.0];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Box,x0Box)

# CoAxCop2
x0Body2 = [[0.16;0.0;0.0];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Body2,x0Body2)

x0Gimbal2 = [[0.16;0.0;0.143];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Gimbal2,x0Gimbal2)

x0Prop21 = [[0.16;0.0;0.234];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop21,x0Prop21)

x0Prop22 = [[0.16;0.0;0.24];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop22,x0Prop22)

## Joint Descriptions
# Joint 1 (Free) between the inertial frame and body.
j1 = Joint(InFrame,Body1,zeros(3),zeros(3))

# CoAxCop1
# Joint 2 (Revolute2) between the body and gimbal
rjBody = [0.0 0.0 0.134][:]
rjGimbal1 = [0.0 0.0 -0.009][:]
j2 = Joint(Body1,Gimbal1,rjBody,rjGimbal1,
     type = "Revolute2")#,jointTorque = [0.0 0.0 0.0][:])

# Joint 3 (Revolute) between the body and gimbal
rjGimbal2 = [0.0 0.0 0.091][:]
rjProp1 = [0.0 0.0 0.0][:]
j3 = Joint(Gimbal1, Prop11, rjGimbal2, rjProp1, type = "Revolute",jointTorque = [0.0 0.0 0.001][:])

# Joint 4 (Revolute) between the body and gimbal
rjGimbal3 = [0.0 0.0 0.0965][:]
rjProp2 = [0.0 0.0 0.0][:]
j4 = Joint(Gimbal1, Prop12, rjGimbal3, rjProp2, type = "Revolute", jointTorque = [0.0 0.0 -0.001][:])

# Between CACs and Box
# Joint 5 between CoAxCop1 and Box
rjCac1 = [0.04 0.0 0.0][:]
rjBox1 = [-0.04 0.0 0.0][:]
j5 = Joint(Body1, Box, rjCac1, rjBox1, type = "Weld")

# Joint 6 between Box and CoAxCop2
rjBox2 = [0.04 0.0 0.0][:]
rjCac2 = [-0.04 0.0 0.0][:]
j6 = Joint(Box, Body2, rjBox2, rjCac2, type = "Weld")

# CoAxCop2
# Joint 7 (Revolute2) between the body and gimbal
rjBody = [0.0 0.0 0.134][:]
rjGimbal1 = [0.0 0.0 -0.009][:]
j7 = Joint(Body2,Gimbal2,rjBody,rjGimbal1,
     type = "Revolute2")#,jointTorque = [0.0 0.0 0.0][:])

# Joint 8 (Revolute) between the body and gimbal
rjGimbal2 = [0.0 0.0 0.091][:]
rjProp1 = [0.0 0.0 0.0][:]
j8 = Joint(Gimbal2, Prop21, rjGimbal2, rjProp1, type = "Revolute",jointTorque = [0.0 0.0 0.001][:])

# Joint 9 (Revolute) between the body and gimbal
rjGimbal3 = [0.0 0.0 0.0965][:]
rjProp2 = [0.0 0.0 0.0][:]
j9 = Joint(Gimbal2, Prop22, rjGimbal3, rjProp2, type = "Revolute", jointTorque = [0.0 0.0 -0.001][:])


## Simulation
tEnd = 0.1
tSpan = 0.01
g = [0.0;0.0;0.0]
j = (j1,j2,j3,j4,j5,j6,j7,j8,j9) # Tuple of joints
tSim, solFinal = simulate(tEnd,tSpan,j...,g=g)

## Assigning solutions
solBody1 = solFinal[1]
solGimbal1 = solFinal[2]
solProp11 = solFinal[3]
solProp12 = solFinal[4]
solBox = solFinal[5]
solBody2 = solFinal[6]
solGimbal2 = solFinal[7]
solProp21 = solFinal[8]
solProp22 = solFinal[9]

## Plotting
jointLoc2 = Matrix{Float64}(undef,length(tSim),3)
jointLoc3 = deepcopy(jointLoc2)
jointLoc4 = deepcopy(jointLoc2)
jointLoc5 = deepcopy(jointLoc2)
jointLoc6 = deepcopy(jointLoc2)
jointLoc7 = deepcopy(jointLoc2)
jointLoc8 = deepcopy(jointLoc2)
jointLoc9 = deepcopy(jointLoc2)

ωBody1 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal1 = Matrix{Float64}(undef,length(tSim),3)
ωProp11 = Matrix{Float64}(undef,length(tSim),3)
ωProp12 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal1InBody1 = Matrix{Float64}(undef,length(tSim),3)
ωRel1 = Matrix{Float64}(undef,length(tSim),3)
ωRel2 = Matrix{Float64}(undef,length(tSim),3)
ωRel3 = Matrix{Float64}(undef,length(tSim),3)

ωBox = Matrix{Float64}(undef,length(tSim),3)

ωBody2 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal2 = Matrix{Float64}(undef,length(tSim),3)
ωProp21 = Matrix{Float64}(undef,length(tSim),3)
ωProp22 = Matrix{Float64}(undef,length(tSim),3)
ωGimbal2InBody2 = Matrix{Float64}(undef,length(tSim),3)
ωRel4 = Matrix{Float64}(undef,length(tSim),3)
ωRel5 = Matrix{Float64}(undef,length(tSim),3)
ωRel6 = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    Body1.dcm = quat2dcm(solBody1.β[i,:])
    Gimbal1.dcm = quat2dcm(solGimbal1.β[i,:])
    Prop11.dcm = quat2dcm(solProp11.β[i,:])
    Prop12.dcm = quat2dcm(solProp12.β[i,:])

    ωBody1[i,:] = angVel(solBody.β[i,:],solBody.βdot[i,:])
    ωGimbal1[i,:] = angVel(solGimbal1.β[i,:],solGimbal1.βdot[i,:])
    ωProp11[i,:] = angVel(solProp11.β[i,:],solProp11.βdot[i,:])
    ωProp12[i,:] = angVel(solProp12.β[i,:],solProp12.βdot[i,:])

    ωGimbal1InBody1[i,:] = quat2dcm(solBody1.β[i,:])*
                         transpose(quat2dcm(solGimbal1.β[i,:]))*
                         ωGimbal1[i,:]
    ωRel1[i,:] = ωGimbal1InBody1[i,:] - ωBody1[i,:]
    ωRel2[i,:] = quat2dcm(solGimbal.β[i,:])*
                 transpose(quat2dcm(solProp1.β[i,:]))*
                 ωProp1[i,:] - ωGimbal[i,:]
    ωRel3[i,:] = quat2dcm(solGimbal.β[i,:])*
                 transpose(quat2dcm(solProp2.β[i,:]))*
                 ωProp2[i,:] - ωGimbal[i,:]

    jointLoc2[i,:] = solBody1.r[i,:] + transpose(Body1.dcm)*rjBody  -                    solGimbal.r[i,:] -                     transpose(Gimbal1.dcm)*rjGimbal1
    jointLoc3[i,:] = solGimbal.r[i,:] + transpose(Gimbal1.dcm)*rjGimbal2 - solProp1.r[i,:] - transpose(Prop11.dcm)*rjProp1
    jointLoc4[i,:] = solGimbal.r[i,:] +
                     transpose(Gimbal1.dcm)*rjGimbal3 -
                     solProp2.r[i,:] - transpose(Prop12.dcm)*rjProp2
end
for i=1:length(tSim)
     Body2.dcm = quat2dcm(solBody2.β[i,:])
     Gimbal2.dcm = quat2dcm(solGimbal2.β[i,:])
     Prop21.dcm = quat2dcm(solProp21.β[i,:])
     Prop22.dcm = quat2dcm(solProp22.β[i,:])

     ωBody2[i,:] = angVel(solBody.β[i,:],solBody.βdot[i,:])
     ωGimbal2[i,:] = angVel(solGimbal2.β[i,:],solGimbal2.βdot[i,:])
     ωProp21[i,:] = angVel(solProp21.β[i,:],solProp21.βdot[i,:])
     ωProp22[i,:] = angVel(solProp22.β[i,:],solProp22.βdot[i,:])

     ωGimbal2InBody2[i,:] = quat2dcm(solBody2.β[i,:])*
                          transpose(quat2dcm(solGimbal2.β[i,:]))*
                          ωGimbal2[i,:]
     ωRel4[i,:] = ωGimbal2InBody2[i,:] - ωBody2[i,:]
     ωRel5[i,:] = quat2dcm(solGimbal.β[i,:])*
                  transpose(quat2dcm(solProp1.β[i,:]))*
                  ωProp1[i,:] - ωGimbal[i,:]
     ωRel6[i,:] = quat2dcm(solGimbal.β[i,:])*
                  transpose(quat2dcm(solProp2.β[i,:]))*
                  ωProp2[i,:] - ωGimbal[i,:]

     jointLoc7[i,:] = solBody2.r[i,:] + transpose(Body2.dcm)*rjBody  -                    solGimbal.r[i,:] -                     transpose(Gimbal2.dcm)*rjGimbal2
     jointLoc8[i,:] = solGimbal.r[i,:] + transpose(Gimbal2.dcm)*rjGimbal2 - solProp1.r[i,:] - transpose(Prop21.dcm)*rjProp1
     jointLoc9[i,:] = solGimbal.r[i,:] +
                      transpose(Gimbal2.dcm)*rjGimbal3 -
                      solProp2.r[i,:] - transpose(Prop22.dcm)*rjProp2
end

## Plotting
# Joint Locations
plotPos(tSim,jointLoc2)
plotPos(tSim,jointLoc3)
plotPos(tSim,jointLoc4)
plotPos(tSim,jointLoc7)
plotPos(tSim,jointLoc8)
plotPos(tSim,jointLoc9)

plotErrNorm(tSim,solBody1.β)
plotErrNorm(tSim,solGimbal1.β)
plotErrNorm(tSim,solProp11.β)
plotErrNorm(tSim,solProp12.β)
plotErrNorm(tSim,solBody2.β)
plotErrNorm(tSim,solGimbal2.β)
plotErrNorm(tSim,solProp21.β)
plotErrNorm(tSim,solProp22.β)

plotAngVel(tSim,ωBody1)
plotAngVel(tSim,ωGimbal1)
plotAngVel(tSim,ωProp11)
plotAngVel(tSim,ωProp12)
plotAngVel(tSim,ωRel1)
plotAngVel(tSim,ωRel2)
plotAngVel(tSim,ωRel3)
plotAngVel(tSim,ωBody2)
plotAngVel(tSim,ωGimbal2)
plotAngVel(tSim,ωProp21)
plotAngVel(tSim,ωProp22)
plotAngVel(tSim,ωRel4)
plotAngVel(tSim,ωRel5)
plotAngVel(tSim,ωRel6)
# plotAngVel(tSim,ωGimbalInBody)
# plotAngVel(tSim,ωCube - ωProp)
