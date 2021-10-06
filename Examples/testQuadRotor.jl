# To test a simplified quadrotor.
include("../src/plotSol.jl")
include("../src/simulate.jl")
# include("../src/linearize.jl")
include("../src/nRotor.jl")
# include("../src/trim_kronLazy.jl");
# include("../src/trim_kron.jl")
# include("../src/trim_FinDiff.jl")


# using nRotor.jl, generate quadrotor joint tree
qRotor_j = gen_nRotor(4);
g =([0.0,0.0,-9.81]); # Gravity Vector.

## Simulate
tEnd = 0.1;
tSpan = 0.01;
tSim, solFinal = simulate(tEnd, tSpan, qRotor_j, g=g);

solQuad = solFinal[1];
solProp1 = solFinal[2];
solProp2 = solFinal[3];
solProp3 = solFinal[4];
solProp4 = solFinal[5];

##
ωCube, ωProp1, ωRel1, jointLoc1 = checkRevJoint(solQuad, solProp1, qRotor_j[2].pos1, qRotor_j[2].pos2);
ωCube, ωProp2, ωRel2, jointLoc2 = checkRevJoint(solQuad, solProp2, qRotor_j[3].pos1, qRotor_j[3].pos2);
ωCube, ωProp3, ωRel3, jointLoc3 = checkRevJoint(solQuad, solProp3, qRotor_j[4].pos1, qRotor_j[4].pos2);
ωCube, ωProp4, ωRel4, jointLoc4 = checkRevJoint(solQuad, solProp4, qRotor_j[5].pos1, qRotor_j[5].pos2);

## plots
# errors
# close("all");
# plotErrNorm(tSim,solQuad.β)
# plotErrNorm(tSim,solProp1.β)
# plotErrNorm(tSim,solProp2.β)
# plotErrNorm(tSim,solProp3.β)
# plotErrNorm(tSim,solProp4.β)
# plotPos(tSim,jointLoc1);
# plotPos(tSim,jointLoc2);
# plotPos(tSim,jointLoc3);
# plotPos(tSim,jointLoc4);
#
# # trajectory
# plotPos(tSim, solQuad.r, "Quad");
# plotPos(tSim, solProp1.r, "Prop1");
# plotPos(tSim, solProp2.r, "Prop2");
# plotPos(tSim, solProp3.r, "Prop3");
# plotPos(tSim, solProp4.r, "Prop4");
#
# # angular velocities
# plotAngVel(tSim, ωCube)
# plotAngVel(tSim, ωProp1)

## find trim point and obtain linear plant
x0, u0 = getXU_0(qRotor_j);
nB = qRotor_j[end].RB2.bodyID;
ix= zeros(Integer,20*(nB-1)+1);
# WITH no fault in actuators, hover condition
ix[14*(nB-1)+1:20*(nB-1)] .= 1;
trim_U = u0;
iy = zeros(Integer,7*(nB-1));

out = trim_kronLazy(qRotor_j,g,ix=ix, iy=iy);
# out = trim_FinDiff(qRotor_j,g,ix=ix, iy=iy);
trim_x, trim_u, gam = out;
println("trim_x = ", trim_x)
trim_X = [qRotor_j[1].RB1.x;trim_x]
# trim_U = [zeros(6) reshape(trim_u,(6,nB-1))];
# println("trim_U = ", trim_U[:,2])
# trim_U = zeros(6,nB);
println("qconstr = ", norm(trim_x[4:7])-1)
println("ẍ =", norm(fxdot(trim_X,trim_U, qRotor_j, g)))

## old
## find linear system
x0, u0 = getXU_0(qRotor_j);
x0_woIn = x0[15:end]; u0_woIn = u0[:,2:end]; # without the inertial frame
println()
println("####Linearizing####")
A, Bu = linearize(x0, u0, qRotor_j, g);
Cy = Matrix{Float64}(I,size(A)...); # full state measurement
Cz = Matrix{Float64}(I,size(A)...); # full state output
Dzu = zeros(size(Bu));
## lti system (nominal, no disturbance or faults)
using PyCall
pyc = pyimport("control");
Usys = zeros(length(u0_woIn),length(tSim));
for i=1:length(tSim)
  Usys[:,i] = u0_woIn[:];
end
sysNom = pyc.ss(A, Bu, Cz, Dzu); # nominal system at hover
T, yout, xout = pyc.forced_response(sysNom, permutedims(tSim), Usys, x0_woIn);

## with disturbance
using Random
rng = MersenneTwister(1234);
ξ = randn(1,length(tSim));
Uξsys = [Usys; ξ];

Dyξ = ones(length(x0_woIn)); # disturbance in all states
# no input in measurement
Bξ = ones(length(x0_woIn)); # disturbance in all states
Buξ = [Bu Bξ];
Dzξ = ones(length(x0_woIn)); # disturbance in all states
Dzuξ = [Dzu Dzξ];
sysDist = pyc.ss(A, Buξ, Cz, Dzuξ); # nominal system at hover
T2, yout2, xout2 = pyc.forced_response(sysDist, permutedims(tSim), Uξsys, x0_woIn);

## with faults and disturbance
flt_act = ones(length(u0_woIn)); # actuator fault vector; 1 indicates no fault
act_pos = 9:6:27; # position of thrust in the u vector
flt_act[act_pos] = [1.0;0.5;1.0;0.5]; # partial fault in rotors 2 and 4
Σa = diagm(flt_act);
# BuΣa = Bu*Σa;
Bfξ = [Bu*Σa Bξ];
Dfξ = [Dzu*Σa Dzξ];

flt_sen = ones(size(Cy,1)); # sensor fault vector
Σs = diagm(flt_sen);
Σs_Cy = Σs*Cy; Σs_Dyξ = Σs*Dyξ;

sysFlt = pyc.ss(A, Bfξ, Cz, Dfξ); # system with fault
Tf, yout_f, xout_f = pyc.forced_response(sysFlt, permutedims(tSim), Uξsys, x0_woIn);

## old
# clearconsole()
# j1 = Nothing
# j2 = Nothing
# j3 = Nothing
# j4 = Nothing
# j5 = Nothing
# ##
# mb = 0.155; mp = 0.02;
#
# Ib =  [ 0.0009    0.0000   -0.0000;
#         0.0000    0.0053    0.0000;
#        -0.0000    0.0000    0.0054]#*(10^-3)
#
# Ip =    [ 0.0029    0.0000    0.0000;
#           0.0000    0.5824    0.0000;
#           0.0000    0.0000    0.5805]*(10^-4)
#
# InFrame = InertialFrameAsRB()
# QuadCube = RigidBody(mb,Ib,2)
# Prop1 = RigidBody(mp,Ip,3)
# Prop2 = RigidBody(mp,Ip,4)
# Prop3 = RigidBody(mp,Ip,5)
# Prop4 = RigidBody(mp,Ip,6)
#
# x0Cube = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
# initialiseRigidBody!(QuadCube,x0Cube)
#
# x0Prop1 = [[0.5;-0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
# initialiseRigidBody!(Prop1,x0Prop1)
# x0Prop2 = [[0.5;0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
# initialiseRigidBody!(Prop2,x0Prop2)
# x0Prop3 = [[-0.5;-0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
# initialiseRigidBody!(Prop3,x0Prop3)
# x0Prop4 = [[-0.5;0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
# initialiseRigidBody!(Prop4,x0Prop4)
#
# axis = [0.0 0.0 1.0][:]
# rjCube1 = [0.5 -0.5 0.1][:]
# rjCube2 = [0.5 0.5 0.1][:]
# rjCube3 = [-0.5 -0.5 0.1][:]
# rjCube4 = [-0.5 0.5 0.1][:]
# rjProp = [0.0 0.0 0.0][:];
#
# j1 = Joint(InFrame,QuadCube,zeros(3),zeros(3));
# j2 = Joint(QuadCube,Prop1,rjCube1,rjProp,type = "Revolute",axis = axis);#,jointTorque = [0.0;0.0;0.0001]);
# j3 = Joint(QuadCube,Prop2,rjCube2,rjProp,type = "Revolute",axis = axis);
# j4 = Joint(QuadCube,Prop3,rjCube3,rjProp,type = "Revolute",axis = axis);
# j5 = Joint(QuadCube,Prop4,rjCube4,rjProp,type = "Revolute",axis = axis);
#
# j = ([j1,j2,j3,j4,j5]);
# ## Linearization
# # x0, u0 = getXU_0(j);
# # A,B = linearizeDiff(x0, u0, j);
# # println("A = ",A)
# # println("B = ",B)
# # sleep(1200)
#
# ## Simulation
# tEnd = 0.1;
# tSpan = 0.01;
# g =([0.0,0.0,0.0]) # Gravity Vector.
# tSim, solFinal = simulateDiff(tEnd,tSpan,j,g=g);
#
# solQuad = solFinal[1];
# solProp1 = solFinal[2];
# solProp2 = solFinal[3];
# solProp3 = solFinal[4];
# solProp4 = solFinal[5];
#
# ##
# ωCube, ωProp1, ωRel1, jointLoc1 = checkRevJoint(solQuad, solProp1, rjCube1, rjProp);
# ωCube, ωProp2, ωRel2, jointLoc2 = checkRevJoint(solQuad, solProp2, rjCube2, rjProp);
# ωCube, ωProp3, ωRel3, jointLoc3 = checkRevJoint(solQuad, solProp3, rjCube3, rjProp);
# ωCube, ωProp4, ωRel4, jointLoc4 = checkRevJoint(solQuad, solProp4, rjCube4, rjProp);
# #
#
# # Check if joint location has moved
# # jointLoc1 = Matrix{Float64}(undef,length(tSim),3);
#
# # ωCube = Matrix{Float64}(undef,length(tSim),3);
# # ωProp1 = Matrix{Float64}(undef,length(tSim),3);
# # ωProp1InCube = Matrix{Float64}(undef,length(tSim),3);
# # ωRel = Matrix{Float64}(undef,length(tSim),3);
# # for i=1:length(tSim)
# #     QuadCube.dcm = quat2dcm(solQuad.β[i,:])
# #     ωCube[i,:] = angVel(solQuad.β[i,:],solQuad.βdot[i,:])
# #     ωProp1[i,:] = angVel(solProp1.β[i,:],solProp1.βdot[i,:])
# #     ωProp1InCube[i,:] = quat2dcm(solQuad.β[i,:])*transpose(quat2dcm(solProp1.β[i,:]))*ωProp1[i,:]
# #     ωRel[i,:] = ωProp1InCube[i,:] - ωCube[i,:]
# #     jointLoc1[i,:] = solQuad.r[i,:] + transpose(QuadCube.dcm)*rjCube1 - solProp1.r[i,:]
# # end
# ##
# close("all");
# # plotErrNorm(tSim,solQuad.β)
# # plotErrNorm(tSim,solProp1.β)
# plotPos(tSim,jointLoc2);
# # plotVel(tSim,solQuad.v  - solProp.v);
# # plotAngVel(tSim,ωCube);
# # plotAngVel(tSim,ωRel1);
# # plotAngVel(tSim,ωProp1)
# # plotQuatDot(tSim, solQuad.βdot)
# ##
# # close("all");
# # plotPos(tSim, solQuad.r, "Quad");
# # # plotQuat(tSim, solQuad.β)
# # # plotQuatDot(tSim, solQuad.βdot)
# # # plotErrNorm(tSim, solQuad.β)
# # plotPos(tSim, solProp1.r, "Prop1");
# # # plotVel(tSim, solProp1.v, "Prop1");
# # plotPos(tSim, solProp2.r, "Prop2");
# # plotPos(tSim, solProp3.r, "Prop3");
# # plotPos(tSim, solProp4.r, "Prop4");
# ##
# plotAngVel(tSim, ωCube)
# plotAngVel(tSim, ωProp1)
# plotAngVel(tSim, ωProp3)
#
