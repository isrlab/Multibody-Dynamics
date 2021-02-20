include("../../../src/plotSol.jl")
include("../../../src/simulate.jl")
include("../../../src/linearize.jl")
include("../../../src/OrientationConversion.jl")
include("../../../src/trim_kronLazy.jl");

using Revise
using BenchmarkTools
## clearing variables
# clearconsole()
j1 = Nothing
j2 = Nothing
##
m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1]

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,2)
R2 = RigidBody(m,I1,3)
RbI = InertialFrameAsRB()

# Perturbed from equilibrium
theta1 = deg2rad(30); theta2 = deg2rad(45);
x1_temp = [l/2*sin(theta1);0.0;-l/2*cos(theta1)];
x2_temp = 2*x1_temp + [l/2*sin(theta2);0.0;-l/2*cos(theta2)];
x0R1 = [x1_temp;[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(R1,x0R1)
x0R2 = [x2_temp;[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(R2,x0R2)

axisY = [0.0 1.0 0.0][:] # Axis about which bar is revolving

rj1 = [0.0 0.0 0.0][:] # Joint Location in body frame of first body
rj2 = -R1.x[1:3]

rj3 = x1_temp; # Joint location in body frame of pend1
rj4 = -x2_temp + 2*x1_temp; # Joint location in body frame of pend2

j1 = Joint(RbI,R1,rj1,rj2,type="Revolute",axis=axisY)
j2 = Joint(R1,R2,rj3,rj4,type="Revolute",axis=axisY)

# External Forces Definition
g = [0.0,0.0,-9.806]; # Gravity Vector.

j = [j1;j2];

## Linearization
# x0, u0 = getXU_0(j);
# A,B = linearize(x0, u0, j, g);

## Simulation
tEnd = 1.0
tSpan = 0.01;

tSim, solFinal = simulateDiff(tEnd,tSpan,j,g=g)

solPend1  = solFinal[1];
solPend2  = solFinal[2];
##
ωPend1, jointLoc1 = checkRevJointIn(solPend1, rj2)
ωPend1, ωPend2, ωRel2, jointLoc2 = checkRevJoint(solPend1, solPend2, rj3, rj4);
##
close("all");
# errors
plotErrNorm(tSim,solPend1.β)
plotErrNorm(tSim,solPend2.β)
plotPos(tSim,jointLoc1);
plotPos(tSim,jointLoc2);

# trajectory
plotPos(tSim, solPend1.r)
plotPos(tSim, solPend2.r)
## Compare with MinReal
r1Min, r2Min, v1Min, v2Min, ω1Min, ω2Min = DoublePendulumMinReal(m,l,theta1,theta2, tInt, tEnd);

# Plotting errors vs MinReal
close("all");
plotPos(tSim,r1Min-solPend1.r)
plotPos(tSim,r2Min-solPend2.r)
plotPos(tSim,v1Min-solPend1.v)
plotPos(tSim,v2Min-solPend2.v)
plotAngVel(tSim,ω1Min-ωPend1)
plotAngVel(tSim,ω2Min-ωPend2)
