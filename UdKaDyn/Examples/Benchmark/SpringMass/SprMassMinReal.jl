include("../../../src/plotSol.jl")
using LinearAlgebra
using DifferentialEquations
using Revise
# using JLD
# clearconsole()
function SprMassMinReal(m,l,k,c,restLen, tInt, tEnd)
    function xdot!(xdot::Vector{Float64},x::Vector{Float64})
        xdot[1] = x[2]
        xdot[2]  = -k/m*(x[1]-restLen) - c/m*x[2];
    end

    function mainDynMinReal!(dx,x,p,t)
        xdot!(dx,x)
    end

    tEnd = 10.0; tSpan = 0.01;
    x0 = [l;1.0]
    prob = ODEProblem(mainDynMinReal!,x0,(0.0,tEnd))
    sol = solve(prob,Tsit5(),saveat = tSpan,reltol=1e-10,abstol=1e-10)
    tSim = sol.t
    xSol = sol[1,:]
    xDotSol = sol[2,:]
    rSol = [xSol+l/2*ones(length(tSim)) zeros(length(tSim),2)]
    vSol = [xDotSol zeros(length(tSim),2)]
    return rSol, vSol
end
