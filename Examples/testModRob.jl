# To test a simplified co ax copter.
# Modeling the system as a cube, gimbal and 2 propellers.
# A revolute2 joint between the cube and gimbal
# 1 Rev Joint each between the gimbal and either propeller.

include("../src/plotSol.jl")
include("../src/OrientationConversion.jl")
include("../src/simulate.jl")
#  include("../src/simulateRK.jl")

using JLD
clearconsole()

## Inertial Properties
const mb = 0.755; const mg = 0.215; const mp = 0.02;

const Ib =  [ 0.0009    0.0000   -0.0000;
        0.0000    0.0053    0.0000;
       -0.0000    0.0000    0.0054]#*(10^-3)

const Ip =    [ 0.0029    0.0000    0.0000;
          0.0000    0.5824    0.0000;
          0.0000    0.0000    0.5805]*(10^-4)

const Ig =    [ 0.1435    0.0000    0.0000;
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
     type = "Revolute",jointTorque = [0.0 0.0 0.0][:])

# Joint 3 (Revolute) between the gimbal and prop11
rjGimbal2 = [0.0 0.0 0.091][:]
rjProp1 = [0.0 0.0 0.0][:]
j3 = Joint(Gimbal1, Prop11, rjGimbal2, rjProp1, type = "Revolute",jointTorque = [0.0 0.0 0.001][:])

# Joint 4 (Revolute) between the gimbal and prop12
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
# const rjBody = [0.0 0.0 0.134][:]
# const rjGimbal1 = [0.0 0.0 -0.009][:]
j7 = Joint(Body2,Gimbal2,rjBody,rjGimbal1,
     type = "Revolute")#,jointTorque = [0.0 0.0 0.0][:])

# Joint 8 (Revolute) between the gimbal and prop21
# const rjGimbal2 = [0.0 0.0 0.091][:]
# const rjProp1 = [0.0 0.0 0.0][:]
j8 = Joint(Gimbal2, Prop21, rjGimbal2, rjProp1, type = "Revolute",jointTorque = [0.0 0.0 0.001][:])

# Joint 9 (Revolute) between the gimbal and prop22
# const rjGimbal3 = [0.0 0.0 0.0965][:]
# const rjProp2 = [0.0 0.0 0.0][:]
j9 = Joint(Gimbal2, Prop22, rjGimbal3, rjProp2, type = "Revolute", jointTorque = [0.0 0.0 -0.001][:])


## Simulation
tEnd = 1.0
tSpan = 0.01
const g = MVector{3}([0.0,0.0,-9.806]) # Gravity Vector.
j = (j1,j2,j3,j4,j5,j6,j7,j8,j9) # Tuple of joints
@time tSim, solFinal = simulate(tEnd,tSpan,j...,g=g);
# @time tSim, solFinal = simulateRK45(tEnd, tSpan, j..., g=g);

# save("ModRob1dSecSimRK45NoTol.jld","solFinal", solFinal,"tSim", tSim)
