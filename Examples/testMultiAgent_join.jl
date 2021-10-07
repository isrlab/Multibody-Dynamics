# testing multiple agents coming together
# 2 agents considered are cubeprop systems
# they attach to a box and fly together
include("../src/plotSol.jl")
include("../src/simulateDiff.jl")
include("../src/OrientationConversion.jl")

clearconsole()
j1 = Nothing
j2 = Nothing
j3 = Nothing
## Agent 1
mb = 0.155; mp = 0.02;

Ib =  [ 0.0009    0.0000   -0.0000;
        0.0000    0.0053    0.0000;
       -0.0000    0.0000    0.0054]#*(10^-3)

Ip =    [ 0.0029    0.0000    0.0000;
          0.0000    0.5824    0.0000;
          0.0000    0.0000    0.5805]*(10^-4)

InFrame = InertialFrameAsRB()
QuadCube_1 = RigidBody(mb,Ib,2)
Props_1 = RigidBody(mp,Ip,3)

x0Cube_1 = [[-0.5;0.0;0.0];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(QuadCube_1,x0Cube_1)
x0Props_1 = [[-0.5;0.0;0.5];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Props_1,x0Props_1)

axis = [0.0 0.0 1.0][:]
rjCube = [0.0 0.0 0.5][:]
rjProp = [0.0 0.0 0.0][:]

j1_1 = Joint(InFrame,QuadCube_1,zeros(3),[0.5;zeros(2)])
j1_2 = Joint(QuadCube_1,Props_1,rjCube,rjProp,type = "Revolute",axis = axis);#,jointTorque = [0.0;0.0;0.0001])

j1 = [j1_1,j1_2];

## Agent 2
QuadCube_2 = RigidBody(mb,Ib,2)
Props_2 = RigidBody(mp,Ip,3)

x0Cube_2 = [[0.5;0.0;0.0];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(QuadCube_2,x0Cube_2)
x0Props_2 = [[0.5;0.0;0.5];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Props_2,x0Props_2)

axis = [0.0 0.0 1.0][:]
rjCube = [0.0 0.0 0.5][:]
rjProp = [0.0 0.0 0.0][:]

j2_1 = Joint(InFrame,QuadCube_2,[-0.5;zeros(2)],zeros(3))
j2_2 = Joint(QuadCube_2,Props_2,rjCube,rjProp,type = "Revolute",axis = axis);#,jointTorque = [0.0;0.0;0.0001])

j2 = [j2_1,j2_2];
## Box
Box = RigidBody(mb,Ib,2);
x0Box = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]; # Box cm is at the inertial frame
initialiseRigidBody!(Box,x0Box);

j3_1 = Joint(InFrame, Box, zeros(3), zeros(3));

j3 = [j3_1];
## Joining bodies together
# join cubeprop with a box
joinBodies!(j1,j3,2,2,[0.25;zeros(2)],[-0.25;zeros(2)]);
# println(length(j1))
println(j1[end].RB2.bodyID)
# join cubepropbox with another cubeprop
joinBodies!(j1,j2,4,2,[0.25;zeros(2)],[-0.25;zeros(2)]);
println(j1[end].RB2.bodyID)

# joining two cubeprops together
# joinBodies!(j1,j2,2,2,[0.5;zeros(2)],[-0.5;zeros(2)]);
# println(j1[3].RB2.bodyID)
# println(j1[5].RB2.bodyID)
