include("../src/plotSol.jl")
using LinearAlgebra
using DifferentialEquations
using Revise
using JLD
clearconsole()
# cd("/home/vish0908/Documents/ModRob/ISRL-MBD")

global m = 1.0
global l = 1.0
global g = 9.806
function xdot(x::Vector{Float64})
    global l, m, g
    thDiff = x[3] - x[1] # theta2 - theta1
    A = zeros(4,4); b = zeros(4)
    A[1,1] = 1; b[1] = x[2]
    A[2,:] = [0 4/3 0 cos(thDiff)/2]; b[2] = 3*g/2/l*sin(x[1]) + sin(thDiff)/2*(x[4]^2)
    A[3,3] = 1; b[3] = x[4]
    A[4,:] = [0 cos(thDiff)/2 0 1/3]; b[4] = g/2/l*sin(x[3]) - sin(thDiff)/2*(x[2]^2)

    xdot = A\b
    return xdot
end

function mainDynMinReal(x,p,t)
    dx = xdot(x)
    println("t = ", t)
    return dx
end

tEnd = 1.0; tSpan = 0.01;
x0 = [pi/2; 0.0; pi/2; 0.0]
prob = ODEProblem(mainDynMinReal,x0,(0.0,tEnd))
sol = solve(prob,Tsit5(),saveat = tSpan,reltol=1e-10,abstol=1e-10)
tSim = sol.t
th1Sol = sol[1,:]
th1DotSol = sol[2,:]
th2Sol = sol[3,:]
th2DotSol = sol[4,:]
x1Sol = zeros(length(tSim)); x2Sol = zeros(length(tSim));
z1Sol = zeros(length(tSim)); z2Sol = zeros(length(tSim));
y1Sol = zeros(length(tSim));
vSol = zeros(length(tSim),3)
r2Sol = zeros(length(tSim),3)
for i=1:length(tSim)
    x1Sol[i] = l/2*sin(th1Sol[i])
    z1Sol[i] = l/2*cos(th1Sol[i])
    vSol[i,:] = [l/2*cos(th1Sol[i])*th1DotSol[i] 0.0 l/2*sin(th1Sol[i])*th1DotSol[i]]
    x2Sol[i] = l*sin(th1Sol[i]) + l/2*sin(th2Sol[i])
    z2Sol[i] = l*cos(th1Sol[i]) + l/2*cos(th2Sol[i])
end
r1Sol = [x1Sol y1Sol z1Sol]
r2Sol = [x2Sol y1Sol z2Sol]
ω1Sol = [y1Sol th1DotSol y1Sol]
ω2Sol = [y1Sol th2DotSol y1Sol]
close("all");
plotAngVel(tSim,ω1Sol)
plotAngVel(tSim,ω2Sol)
plotPos(tSim, r1Sol)
plotPos(tSim, r2Sol)
##
save("DoublePendulumMin.jld","r1Min",r1Sol,"r2Min", r2Sol)#,"vMin",vSol)

# r1Sol = transpose(sol[1:3,:])
# vSol = transpose(sol[4:6,:])
# eaSol = transpose(sol[7:9,:])
# ωSol = transpose(sol[10:12,:])
#
#
