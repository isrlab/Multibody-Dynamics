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
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end
    @show X0
    # Declaring the ODE Problem as per DifferentialEquations convention
    prob = ODEProblem(mainDynODE!,X0,(0.0,tEnd),j)
    sol = solve(prob,saveat=tSpan,Tsit5())
end

function mainDynODE!(dX,X,j::Tuple{Joint},t)
    # ODE function to be used as per DifferentialEquations covnention
    # Create extForcesList storing extForces for each rigid body
    # Create ForceConstraints Array storing constraint forces acting on each rigid body
    GravityInInertial = [0.0;0.0;-9.806]
    extFList = Vector{extForces}(undef,length(j)+1)
    ForceConstr = Array{Float64}(undef,7,length(j)+1)
    extFList[1] = getExternalForce(j[1].RB1)
    for k=2:length(j)+1
        extFList[k] = getExternalForce(j[k-1].RB2)
    end
    ForceConstr[:,1] = Vector{Float64}(undef,7)
    # ForceConstr[:,1] = ForceCon(j[1],extFList[1],extFList[2],GravityInInertial)
    for k=2:length(j)+1
        ForceConstr[:,k] = ForceCon(j[k-1],extFList[k-1],extFList[k],GravityInInertial)
    end
    mainDyn!(dX,X,j,extFList,ForceConstr)
    # extF1 = extForces(zeros(1,3),zeros(1,3),zeros(1,3))
    # extF2 = extForces(zeros(1,3),zeros(1,3),zeros(1,3))
end

function mainDyn!(dQ,Q,j::Tuple{Joint},extFList::Vector{extForces}, ForceConstr::Array{Float64,2})
    GravityInInertial = [0.0;0.0;-9.806]
    @show Q
    j[1].RB1.x = Q[1:14]
    for k=1:length(j)
        for m=1:14
        j[k].RB2.x[m] = Q[14*k + m]
        end
    end
    # First body always the inertial frame
    dQ[1:14] = zeros(14)
    # rbDynQuat!(dQ[1:14],j[1].RB1,extFList[1],ForceConstr[:,1])
    for k=1:length(j)
        rbDynQuat!(dQ[14*k+1:14*(k+1)],j[k].RB2,extFList[k+1],ForceConstr[:,k+1],GravityInInertial)
    end

    # j.RB1.x = Q[1:14]
    # j.RB2.x = Q[15:28]
    # # Force generated from Constraint (from UK formulation)
    # Fc = ForceCon(j,extF1,extF2,GravityInInertial)
    # Fc1 = Fc[1:7]
    # Fc2 = Fc[8:end]
    # x1dot = dQ[1:14]
    # x2dot = dQ[15:28]
    # rbDynQuat(x1dot,j.RB1,extF1,Fc1,GravityInInertial)
    # rbDynQuat(x2dot,j.RB2,extF2,Fc2,GravityInInertial)
end

function getExternalForce(b::RigidBody)
    # External forces function to be changed later.
    Forces = zeros(1,3)
    Positions = zeros(1,3)
    Torques = zeros(1,3)
    return extF = extForces(Forces,Positions,Torques)
end


# function f(x::Float64...)
#     println(x[1])
#     println(length(x))
#     return sum(x)
# end
#
# f(1.0,2.0,3.0)
