# simulate.jl

```@docs
simulate(tEnd::Float64,tSpan::Float64,j::Joint...;g::Vector{Float64}=[0.0;0.0;-9.806],extFVec::Vector{extForces}=Vector{extForces}(undef,1))

mainDynODE(X::Vector{Float64},j::Tuple{Vararg{Joint}},t::Float64)

mainDyn(Q::Vector{Float64},j::Tuple{Vararg{Joint}},extFList::Vector{extForces}, ForceConstr::Array{Float64,2}, GravityInInertial::Vector{Float64})
```
