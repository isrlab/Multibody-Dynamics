# Native implementation of ode45 to check for speed.
include("Joint.jl")

function rk4(tEnd::Float64, tStart::Float64, x0::Vector{Float64}, p::Tuple{Vararg{Joint}}, f::Function; h::Float64 = 0.01)
    # x0's type has not been declared because we do not know if it is a vector or a single float object.
    # eps: Tolerance, default value: 1e-3
    # To solve xdot = f(t,x)
    # p: Parameters
    # h: tSpan (interval time)

    tSim = Float64[]; push!(tSim,tStart)
    xSol = Float64[]; xSol = vcat(xSol,x0)
    len = length(x0) # Length of vector
    x = x0;
    t = tStart
    while t < tEnd
        k1 = h*f(t, x, p)
        k2 = h*f(t + h/2, x + k1/2, p)
        k3 = h*f(t + h/2, x + k2/2, p)
        k4 = h*f(t + h, x + k3, p)

        y1 = x + 1/6*(k1 + 2*k2 + 2*k3 + k4)

        t = t + h
        push!(tSim,t)
        x = y1
        xSol = vcat(xSol,x)
    end
    return tSim, xSol
end
