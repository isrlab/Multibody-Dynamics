include("simulate.jl")
include("plotSol.jl")
include("OrientationConversion.jl")
using Revise
using JLD
clearconsole()

m = 1.0; l = 1.0; # Mass and length of bar
# Assuming bar revolves about Y axis
I1 = [1 0 0; 0 m*l^2/12 0; 0 0 1]

# Testing Dynamics with Revolute Joint
R1 = RigidBody(m,I1,"quaternions")
RbI = InertialFrameAsRB()
x0R1 = [l/2;zeros(2);[1;zeros(3)];zeros(3);zeros(4)]
R1 = initialiseRigidBody(R1,x0R1)

axis = [0.0 1.0 0.0][:] # Axis about which bar is revolving

rj1 = [0.0 0.0 0.0][:] # Joint Location in body frame of first body
rj2 = [-0.5 0.0 0.0][:] # Joint location in body frame of second body

j = Joint(RbI,R1,rj1,rj2,"Revolute",axis)

# Simulation
tEnd = 1.0
tSpan = 0.01
sol = simulate(tEnd,tSpan,j)

tSim = sol.t
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

# Plotting errors
plotPos(tSim,rMin-rSol)
plotVel(tSim,vMin-vSol)
plotAngVel(tSim,ωMin+ωSol)
