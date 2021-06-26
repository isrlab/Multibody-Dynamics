# include("../src/plotSol.jl")
using LinearAlgebra
using DifferentialEquations
using Revise
# clearconsole()
function PendulumMinReal(m,l,theta,tInt,tEnd)
    function thdot(th)
        thDot = zeros(size(th));
        thDot[1] = th[2]
        thDot[2]  = -3*9.806/2/l*sin(th[1])
        return thDot
    end

    function mainDynMinReal(x,p,t)
        dx = thdot(x);
        return dx
    end

    x0 = [theta; 0.0];
    prob = ODEProblem(mainDynMinReal,x0,(0.0,tEnd))
    sol = solve(prob,Tsit5(),saveat = tInt,reltol=1e-10,abstol=1e-10);
    tSim = sol.t;
    thSol = sol[1,:];
    thDotSol = sol[2,:];
    xSol = zeros(length(tSim))
    zSol = zeros(length(tSim))
    ySol = zeros(length(tSim))
    vSol = zeros(length(tSim),3)
    for i=1:length(tSim)
        xSol[i] = l/2*sin(thSol[i])
        zSol[i] = -l/2*cos(thSol[i])
        vSol[i,:] = [l/2*cos(thSol[i])*thDotSol[i] 0.0 l/2*sin(thSol[i])*thDotSol[i]]
    end
    rSol = [xSol ySol zSol];
    ωSol = [ySol -thDotSol ySol];
    return rSol, vSol, ωSol
end