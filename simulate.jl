# Simulation
# Tentatively to simulate two rigid bodies constrained by one revolute joint


include("RigidBody.jl")
include("OrientationConversion.jl")
include("Joint.jl")
include("Force.jl")


using LinearAlgebra
using Revise
using DifferentialEquations

function simulate(tEnd,tSpan,j::Joint...)
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function

    # Initial condition (all bodies connected)
    X0 = j[1].RB1.x
    for each k= 1:length(j)
        push!(X0,j[k].RB2.x)
    end
    @show X0
    # Declaring the ODE Problem as per DifferentialEquations convention
    prob = ODEProblem(mainDynODE!,X0,(0.0,tEnd),j)
    sol = solve(prob,saveat=tSpan,Tsit5())
end

function mainDynODE!(dX,X,j::Tuple{Joint},t)
    # ODE function to be used as per DifferentialEquations covnention
    # Create extForcesList storing extForces for each rigid body
    extFList = Vector{extForces}(undef,length(j)+1)
    ForceConstr = Array{Float64}(undef,7,length(j)+1)
    for k=1:length(j)+1
        push!(extFList,genExtF(j[k].RB2))

    end

    mainDyn!(dX,X,j,extFList)

    # extF1 = extForces(zeros(1,3),zeros(1,3),zeros(1,3))
    # extF2 = extForces(zeros(1,3),zeros(1,3),zeros(1,3))

end

function mainDyn!(dQ,Q,j::Tuple{Joint},extFList::Vector{extForces})
    j.RB1.x = Q[1:14]
    j.RB2.x = Q[15:28]
    GravityInInertial = [0.0;0.0;-9.806]
    # Force generated from Constraint (from UK formulation)
    Fc = ForceCon(j,extF1,extF2,GravityInInertial)
    Fc1 = Fc[1:7]
    Fc2 = Fc[8:end]
    x1dot = dQ[1:14]
    x2dot = dQ[15:28]
    rbDynQuat(x1dot,j.RB1,extF1,Fc1,GravityInInertial)
    rbDynQuat(x2dot,j.RB2,extF2,Fc2,GravityInInertial)

end


function genMatM(b::RigidBody)
    # Function to generate mass matrix for each rigid body
    β = b.x[4:7]
    E = genE(β)
    J = b.J
    M = [b.m*Matrix{Float64}(I,3,3)          zeros(3,4)
                         zeros(4,3) 4*transpose(E)*J*E]
end

function genExtF(b::RigidBody,extF::extForces,GravityInInertial::Vector{Float64})
    # Function to generate external Force vector for unconstrained system
    β = b.x[4:7]
    βdot = b.x[11:14]
    E = genE(β)
    Edot = -genE(βdot)
    TotalMoment = zeros(3)
    J0 = # random positive number
    for i in 1:size(extF.Forces)[1]
        TotalMoment = TotalMoment + cross(extF.Positions[i,:],extF.Forces[i,:])
    end
    TotalMoment = TotalMoment + sum(extF.Torques,dims=1)[:]
    TotalMoment4 = [0.0;TotalMoment]
    F = [transpose(b.dcm)*(sum(extF.Forces,dims=1)[:]) + b.m*GravityInInertial
         TotalMoment4 - 8*transpose(Edot)*b.J*E*βdot - 4*b.J[1,1]*(transpose(βdot)*βdot)*β]
    return F
end

# function f(x::Float64...)
#     println(x[1])
#     println(length(x))
#     return sum(x)
# end
#
# f(1.0,2.0,3.0)
