# Native implementation of ode45 to check for speed.
include("Joint.jl")

function rk45(tEnd::Float64, tStart::Float64, x0::Vector{Float64}, p::Tuple{Vararg{Joint}}, f::Function; eps::Float64=1e-3)
    # eps: Tolerance, default value: 1e-3
    # To solve xdot = f(t,x)
    # p: Parameters

    tSim = Float64[]; push!(tSim,tStart)
    xSol = Float64[]; xSol = vcat(xSol,x0)
    h = 0.01 # Initial step size. Arbitrary
    len = length(x0) # Length of vector
    x = x0;
    t = tStart
    while t < tEnd
        k1 = h*f(t, x, p)
        k2 = h*f(t + h/4, x + k1/4, p)
        k3 = h*f(t + 3*h/8, x + 3/32*k1 + 9/32*k2, p)
        k4 = h*f(t + 12*h/13, x + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3, p)
        k5 = h*f(t + h, x + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4, p)
        k6 = h*f(t + h/2, x - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5, p)

        y1 = x + 25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - k5/5
        y2 = x + 16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6

        err = y1 - y2
        R = 1/h*broadcast(abs,err)
        δ = 0.84*(eps/norm(R,Inf))^0.25

        if R <= eps*ones(len)
            t = t + h
            x = y1
            xSol = vcat(xSol,x)
            push!(tSim,t)
            h = δ*h
        else
            h = δ*h
        end
    end
    return tSim, xSol
end
