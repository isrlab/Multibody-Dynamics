# Simulation
# Tentatively to simulate two rigid bodies constrained by one revolute joint


include("RigidBody.jl")
include("OrientationConversion.jl")
include("Joint.jl")
include("Force.jl")


using LinearAlgebra
using Revise
using DifferentialEquations

global GravityInInertial = Vector{Float64}(undef,3)

struct solSim
    # Structure to store results of Simulation for each rigid body
    # t::Vector{Float64}
    r::Matrix{Float64}
    v::Matrix{Float64}
    β::Matrix{Float64}
    βdot::Matrix{Float64}
end

function simulate(tEnd::Float64,tSpan::Float64,j::Joint...;g::Vector{Float64}=[0.0;0.0;-9.806])
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function
    # Initial condition (all bodies connected)
    global GravityInInertial = g
    X0 = j[1].RB1.x
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end
    # Declaring the ODE Problem as per DifferentialEquations convention
    prob = ODEProblem(mainDynODE,X0,(0.0,tEnd),j)
    sol = solve(prob,saveat=tSpan,Tsit5(),reltol=1e-10,abstol=1e-10)
    # return sol

    solFinal = Vector{solSim}(undef,length(j))
    for i=1:length(j)
        rSol = transpose(sol[14*i+1:14*i+3,:])
        vSol = transpose(sol[14*i+8:14*i+10,:])
        βSol = transpose(sol[14*i+4:14*i+7,:])
        βdotSol = transpose(sol[14*i+11:14*i+14,:])
        solFinal[i] = solSim(rSol,vSol,βSol,βdotSol)
    end
    return (sol.t, solFinal)
end

function mainDynODE(X::Vector{Float64},j::Tuple{Vararg{Joint}},t::Float64)
    # ODE function to be used as per DifferentialEquations covnention
    # Create extForcesList storing extForces for each rigid body
    # Create ForceConstraints Array storing constraint forces acting on each rigid body
    global GravityInInertial

    # Update RigidBodies
    updateRigidBody(j[1].RB1,X[1:14])
    for k=1:length(j)
        updateRigidBody(j[k].RB2,X[14*k+1:14*(k+1)])
    end

    # Generate extForces for each body
    extFList = Vector{extForces}(undef,length(j)+1)
    extFList[1] = getExternalForce(j[1].RB1)
    for k=2:length(j)+1
        extFList[k] = getExternalForce(j[k-1].RB2)
    end

    # Generate constraint Forces for each body using UK Formulation
    ForceConstr = Array{Float64}(undef,7,length(j)+1)
    ForceConstr[:,2] = ForceCon(j[1],extFList[1],extFList[2],GravityInInertial)
    if length(j) > 1
        for k=2:length(j)
            ForceConstr[:,k] = ForceCon(j[k],extFList[k],extFList[k+1],GravityInInertial)[1:7]
            ForceConstr[:,k+1] = ForceCon(j[k],extFList[k],extFList[k+1],GravityInInertial)[8:14]
        end
    end


    dX = mainDyn(X,j,extFList,ForceConstr, GravityInInertial)
    return dX
end

function mainDyn(Q::Vector{Float64},j::Tuple{Vararg{Joint}},extFList::Vector{extForces}, ForceConstr::Array{Float64,2}, GravityInInertial::Vector{Float64})
    dQ = Vector{Float64}(undef,(length(j)+1)*14)

    # First body always the inertial frame
    dQ[1:14] = zeros(14)
    for k=1:length(j)
        # rbDynQuat!(dQ[14*k+1:14*(k+1)],j[k].RB2,extFList[k+1],ForceConstr[:,k+1],GravityInInertial)
        dQ[14*k+1:14*(k+1)] = rbDynQuat(j[k].RB2,extFList[k+1],ForceConstr[:,k+1],GravityInInertial)
    end
    return dQ
end

function getExternalForce(b::RigidBody)
    # External forces function to be changed later.
    Forces = zeros(1,3)
    Positions = zeros(1,3)
    Torques = zeros(1,3)
    return extF = extForces(Forces,Positions,Torques)
end
