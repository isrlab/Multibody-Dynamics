include("simulate.jl")
include("plotSol.jl")
include("OrientationConversion.jl")

clearconsole()
m = 1.0; l = 1.0;
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1]
using Revise
using JLD
# using Plots
# pyplot();
# Test the code.
# for i in 1:100
#     #ang0 = [45*pi/180;60*pi/180;30*pi/180];
#     ang0 = rand(3);
#     q = angle2quat(ang0,"321");
#     ang = quat2angle(q,"321");
#     e = ang - ang0;
#     println("Norm error:", norm(e));
# end

# Simulate rigid body dynamics
# We have to now call ode solver with parameters
# -- this will have a function that will be called to compute all time varying
# -- quantities.


# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,"quaternions")
RbI = InertialFrameAsRB()
x0R1 = [0.5;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
R1 = initialiseRigidBody(R1,x0R1)
axis = [0.0 1.0 0.0][:]
rj1 = [0.0 0.0 0.0][:]
rj2 = [-0.5 0.0 0.0][:]
# x0R2 = rand(14)
# R2 = RigidBody(m,I1,"quaternions")
# R2 = initialiseRigidBody(R1,x0R2)
j = Joint(RbI,R1,rj1,rj2,"Revolute",axis)
tEnd = 1.0
tSpan = 0.01
sol = simulate(tEnd,tSpan,j)

tSim = LinRange(0.0,tEnd,Int64(tEnd/tSpan)+1)
rSol = transpose(sol[15:17,:])
vSol = transpose(sol[22:24,:])
βSol = transpose(sol[18:21,:])
βdotSol = transpose(sol[25:28,:]);

# Check errNorm
plotErrNorm(tSim,βSol)
# Check if joint location has moved
jointLoc = Matrix{Float64}(undef,length(tSim),3)
ωSol = Matrix{Float64}(undef,length(tSim),3)
for i=1:length(tSim)
    R1.dcm = quat2dcm(βSol[i,:])
    ωSol[i,:] = angVel(βSol[i,:],βdotSol[i,:])
    jointLoc[i,:] = rSol[i,:] + transpose(R1.dcm)*rj2
end
plotPos(tSim,jointLoc)

# Compare with MinReal
solMin = load("PendulumMin.jld")
ωMin = solMin["ωMin"]
rMin = solMin["rMin"]
vMin = solMin["vMin"]
plotPos(tSim,rSol-rMin)
plotVel(tSim,vSol-vMin)
plotAngVel(tSim,ωSol+ωMin)
