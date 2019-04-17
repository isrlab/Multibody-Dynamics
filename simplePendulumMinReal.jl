using LinearAlgebra
using DifferentialEquations
using Revise
using JLD
clearconsole()
# cd("/home/vish0908/Documents/ModRob/ISRL-MBD")

global m = 1.0
global l = 1.0
function thdot(th::Vector{Float64})
    global l, m
    thdot = zeros(2)
    thdot[1] = th[2]
    thdot[2]  = -3*9.806/2/l*sin(th[1])
    return thdot
end

function mainDynMinReal(x,p,t)
    dx = thdot(x)
    return dx
end

tEnd = 1.0; tSpan = 0.01;
x0 = [pi/2;0.0]
prob = ODEProblem(mainDynMinReal,x0,(0.0,tEnd))
sol = solve(prob,Tsit5(),saveat = tSpan,reltol=1e-10,abstol=1e-10)
tSim = sol.t
thSol = sol[1,:]
thDotSol = sol[2,:]
xSol = zeros(length(tSim))
zSol = zeros(length(tSim))
ySol = zeros(length(tSim))
for i=1:length(tSim)
    xSol[i] = l/2*sin(thSol[i])
    zSol[i] = -l/2*cos(thSol[i])
end
rSol = [xSol ySol zSol]
ωSol = [ySol thDotSol ySol]
save("PendulumMin.jld","rMin",rSol,"ωMin", ωSol)

# rSol = transpose(sol[1:3,:])
# vSol = transpose(sol[4:6,:])
# eaSol = transpose(sol[7:9,:])
# ωSol = transpose(sol[10:12,:])
#
#