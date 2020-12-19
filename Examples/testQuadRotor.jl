# To test a simplified quadrotor.

include("../src/plotSol.jl")
include("../src/simulate.jl")
include("../src/OrientationConversion.jl")
##
clearconsole()
j1 = Nothing
j2 = Nothing
j3 = Nothing
j4 = Nothing
j5 = Nothing

##
mb = 0.155; mp = 0.02;

Ib =  [ 0.0009    0.0000   -0.0000;
        0.0000    0.0053    0.0000;
       -0.0000    0.0000    0.0054]#*(10^-3)

Ip =    [ 0.0029    0.0000    0.0000;
          0.0000    0.5824    0.0000;
          0.0000    0.0000    0.5805]*(10^-4)

InFrame = InertialFrameAsRB()
QuadCube = RigidBody(mb,Ib,2)
Prop1 = RigidBody(mp,Ip,3)
Prop2 = RigidBody(mp,Ip,4)
Prop3 = RigidBody(mp,Ip,5)
Prop4 = RigidBody(mp,Ip,6)

x0Cube = [zeros(3);[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(QuadCube,x0Cube)

x0Prop1 = [[0.5;-0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop1,x0Prop1)
x0Prop2 = [[0.5;0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop2,x0Prop2)
x0Prop3 = [[-0.5;-0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop3,x0Prop3)
x0Prop4 = [[-0.5;0.5;0.1];[1;zeros(3)];zeros(3);zeros(4)]
initialiseRigidBody!(Prop4,x0Prop4)

axis = [0.0 0.0 1.0][:]
rjCube1 = [0.5 -0.5 0.1][:]
rjCube2 = [0.5 0.5 0.1][:]
rjCube3 = [-0.5 -0.5 0.1][:]
rjCube4 = [-0.5 0.5 0.1][:]
rjProp = [0.0 0.0 0.0][:];

j1 = Joint(InFrame,QuadCube,zeros(3),zeros(3));
j2 = Joint(QuadCube,Prop1,rjCube1,rjProp,type = "Revolute",axis = axis);#,jointTorque = [0.0;0.0;0.0001]);
j3 = Joint(QuadCube,Prop2,rjCube2,rjProp,type = "Revolute",axis = axis);
j4 = Joint(QuadCube,Prop3,rjCube3,rjProp,type = "Revolute",axis = axis);
j5 = Joint(QuadCube,Prop4,rjCube4,rjProp,type = "Revolute",axis = axis);

## Simulation
tEnd = 0.1;
tSpan = 0.01;
g = MVector{3}([0.0,0.0,0.0]) # Gravity Vector.
tSim, solFinal = simulate(tEnd,tSpan,j1,j2,j3,j4,j5,g=g);

solQuad = solFinal[1];
solProp1 = solFinal[2];
solProp2 = solFinal[3];
solProp3 = solFinal[4];
solProp4 = solFinal[5];

##
ωCube, ωProp1, ωRel1, jointLoc1 = checkRevJoint(solQuad, solProp1, rjCube1, rjProp);
ωCube, ωProp2, ωRel2, jointLoc2 = checkRevJoint(solQuad, solProp2, rjCube2, rjProp);
ωCube, ωProp3, ωRel3, jointLoc3 = checkRevJoint(solQuad, solProp3, rjCube3, rjProp);
ωCube, ωProp4, ωRel4, jointLoc4 = checkRevJoint(solQuad, solProp4, rjCube4, rjProp);
#

# Check if joint location has moved
# jointLoc1 = Matrix{Float64}(undef,length(tSim),3);

# ωCube = Matrix{Float64}(undef,length(tSim),3);
# ωProp1 = Matrix{Float64}(undef,length(tSim),3);
# ωProp1InCube = Matrix{Float64}(undef,length(tSim),3);
# ωRel = Matrix{Float64}(undef,length(tSim),3);
# for i=1:length(tSim)
#     QuadCube.dcm = quat2dcm(solQuad.β[i,:])
#     ωCube[i,:] = angVel(solQuad.β[i,:],solQuad.βdot[i,:])
#     ωProp1[i,:] = angVel(solProp1.β[i,:],solProp1.βdot[i,:])
#     ωProp1InCube[i,:] = quat2dcm(solQuad.β[i,:])*transpose(quat2dcm(solProp1.β[i,:]))*ωProp1[i,:]
#     ωRel[i,:] = ωProp1InCube[i,:] - ωCube[i,:]
#     jointLoc1[i,:] = solQuad.r[i,:] + transpose(QuadCube.dcm)*rjCube1 - solProp1.r[i,:]
# end
##
close("all");
# plotErrNorm(tSim,solQuad.β)
# plotErrNorm(tSim,solProp1.β)
plotPos(tSim,jointLoc2);
# plotVel(tSim,solQuad.v  - solProp.v);
# plotAngVel(tSim,ωCube);
# plotAngVel(tSim,ωRel1);
# plotAngVel(tSim,ωProp1)
# plotQuatDot(tSim, solQuad.βdot)
##
# close("all");
plotPos(tSim, solQuad.r, "Quad");
# plotQuat(tSim, solQuad.β)
# plotQuatDot(tSim, solQuad.βdot)
# plotErrNorm(tSim, solQuad.β)
plotPos(tSim, solProp1.r, "Prop1");
# plotVel(tSim, solProp1.v, "Prop1");
plotPos(tSim, solProp2.r, "Prop2");
plotPos(tSim, solProp3.r, "Prop3");
plotPos(tSim, solProp4.r, "Prop4");


# using JLD
# save("Quad50s.jld","tSim",tSim,"solFinal", solFinal)


# plotAngVel(tSim,ωRel)
# plotAngVel(tSim,ωCube - ωProp)

# # Check revJoint axis
# axisSol = Matrix{Float64}(undef,length(tSims),3)
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
