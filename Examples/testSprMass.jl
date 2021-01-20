include("../src/plotSol.jl")
include("../src/simulateDiff.jl")
include("../src/OrientationConversion.jl")
using Revise
using JLD
clearconsole()

# To test a Spring Mass system
# Impulse of 1.0 m/s in x, no gravity
m = 1.0; l = 1.0; # Mass and length of bar
I1 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,2)
RbI = InertialFrameAsRB()
x0R1 = [1.5*l;zeros(2);[1;zeros(3)];[1;zeros(2)];zeros(4)]
initialiseRigidBody!(R1,x0R1)

rj1 = [0.0 0.0 0.0][:] # Joint Location in body frame of first body
rj2 = [-l/2 0.0 0.0][:] # Joint location in body frame of second body
k = 10.0
j1 = Joint(RbI,R1,rj1,rj2,type="Spring",k=k, rL = l);
j = [j1];

# Simulation
tEnd = 10.0
tSpan = 0.01
g = [0.0;0.0;0.0];
tSim, solFinal = simulateDiff(tEnd,tSpan,j,g = g)

# Check errNorm of quat
plotErrNorm(tSim,solFinal[1].β)

# Compare with MinReal
solMin = load("SprMassMin.jld")
rMin = solMin["rMin"]
vMin = solMin["vMin"]

# Plotting errors
plotPos(tSim,rMin-solFinal[1].r)
plotVel(tSim,vMin-solFinal[1].v)

# tSim = sol.t
# rSol = transpose(sol[15:17,:])
# vSol = transpose(sol[22:24,:])
# βSol = transpose(sol[18:21,:])
# βdotSol = transpose(sol[25:28,:]);
