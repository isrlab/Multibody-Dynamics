# Simulation using native rk45 implementation
# Tentatively to simulate two rigid bodies constrained by one revolute joint


include("RigidBody.jl")
include("OrientationConversion.jl")
include("Joint.jl")
include("Force.jl")
include("extF.jl")
include("rk45.jl")
include("rk4.jl")

using LinearAlgebra
using Revise
using DifferentialEquations
using StaticArrays

global GravityInInertial = MVector{3}([0.0,0.0,-9.806])
# global extFList = Vector{extForces}(undef,2)

struct solSim
    # Structure to store results of Simulation for each rigid body
    r::Matrix{Float64}
    v::Matrix{Float64}
    β::Matrix{Float64}
    βdot::Matrix{Float64}
end

## RK45
function simulateRK45(tEnd::Float64,tSpan::Float64,j::Joint...;
                  g::MArray{Tuple{3},Float64,1,3}=[0.0;0.0;-9.806])#,extFVec::Vector{extForces}=Vector{extForces}(undef,1))
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function
    # Initial condition (all bodies connected)
    global GravityInInertial = g
    # global extFList = extFVec
    X0 = j[1].RB1.x
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end

    # Using native RK implementation
    tSim, sol = rk45(tEnd,0.0,X0,j,mainDynRK,eps = 1e-6)
    nBodies = length(j) + 1
    ntStops = Int64(length(sol)/(14*nBodies))
    sol = reshape(sol,(14*nBodies,ntStops))
    @show sol
    error("break.")
    solFinal = Vector{solSim}(undef,length(j))
    rSol = Matrix{Float64}(undef,ntStops,3)
    vSol = Matrix{Float64}(undef,ntStops,3)
    βSol = Matrix{Float64}(undef,ntStops,4)
    βdotSol = Matrix{Float64}(undef,ntStops,4)
    for i=1:length(j)
        for t=1: ntStops
            rSol[t,:] = sol[14*i+1:14*i+3,t]
            vSol[t,:] = sol[14*i+8:14*i+10,t]
            βSol[t,:] = sol[14*i+4:14*i+7,t]
            βdotSol[t,:] = sol[14*i+11:14*i+14,t]
        end
        solFinal[i] = solSim(rSol,vSol,βSol,βdotSol)
    end
    return tSim, solFinal
end

## RK4
function simulateRK4(tEnd::Float64,tSpan::Float64,j::Joint...;
                  g::MArray{Tuple{3},Float64,1,3}=[0.0;0.0;-9.806])#,extFVec::Vector{extForces}=Vector{extForces}(undef,1))
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function
    # Initial condition (all bodies connected)
    global GravityInInertial = g
    # global extFList = extFVec
    X0 = j[1].RB1.x
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end

    # Using native RK implementation
    sol = rk4(tEnd,0.0,X0,j,mainDynRK, h = tSpan)
    nBodies = length(j) + 1
    ntStops = Int64(length(sol)/(14*nBodies))
    tSim, sol = reshape(sol,(14*nBodies,ntStops))

    solFinal = Vector{solSim}(undef,length(j))
    rSol = Matrix{Float64}(undef,ntStops,3)
    vSol = Matrix{Float64}(undef,ntStops,3)
    βSol = Matrix{Float64}(undef,ntStops,3)
    βdotSol = Matrix{Float64}(undef,ntStops,3)
    for t=1: ntStops
        for i=1:length(j)
            rSol[t,:] = sol[14*i+1:14*i+3,t]
            vSol[t,:] = sol[14*i+8:14*i+10,t]
            βSol[t,:] = sol[14*i+4:14*i+7,t]
            βdotSol[t,:] = sol[14*i+11:14*i+14,t]
        end
        solFinal[i] = solSim(rSol,vSol,βSol,βdotSol)
    end
    return tSim, solFinal
end

function mainDynRK(t::Float64, x::Vector{Float64}, p::Tuple{Vararg{Joint}})
    dx = Vector{Float64}(undef,length(x))
    return mainDynODE!(dx,x,p,t)
end

function mainDynODE!(dX::Vector{Float64}, X::Vector{Float64}, j::Tuple{Vararg{Joint}}, t::Float64)
    # ODE function to be used as per DifferentialEquations convention
    # Create extForcesList storing extForces for each rigid body
    # Create ForceConstraints Array storing constraint forces acting on each rigid body
    global GravityInInertial
    # Update RigidBodies
    updateRigidBody!(j[1].RB1,X[1:14])
    for k=1:length(j)
        updateRigidBody!(j[k].RB2,X[14*k+1:14*(k+1)])
    end

    extFListCopy = externalFTotal(t,j) # Generates all the external forces explicitly specified, through joint actions or directly.

    unconstrF, constrF = Constraint(j, extFListCopy, GravityInInertial)
    # unconstrF is the matrix of the unconstrained Forces acting on the robot.
    # constrF is the matrix of the constraint Forces acting on the system developed using UK formulation

    dX[1:14] = zeros(14)
    for k=1:length(j)
        dX[14*k+1:14*(k+1)] = rbDynQuat(j[k].RB2, unconstrF[:,j[k].RB2.bodyID], constrF[:,j[k].RB2.bodyID])
    end
    # mainDyn!(dX, X, j, unconstrF, constrF)
    return dX
end

function externalFTotal(t::Float64, j::Tuple{Vararg{Joint}})
    # Returns externally applied forces in total, not including gravity. Includes Forces from joints and other explicitly specified external forces.
    extFList = extF(t,j)

    # Generate forces from actuated joints on each body
    ForceJoints = Matrix{Float64}(undef,6,2*length(j))
    for k=1:length(j)
        ForceJoints[:,2*k-1], ForceJoints[:,2*k] = genJointF(t,j[k])
    end

    # Add ForceJoints to extFList
    extFListCopy = deepcopy(extFList)
    for k=1:length(j)
        extFListCopy[j[k].RB1.bodyID].Forces =
        vcat(extFListCopy[j[k].RB1.bodyID].Forces, reshape(ForceJoints[1:3,2*k-1],(1,3)))
        extFListCopy[j[k].RB1.bodyID].Positions =
        vcat(extFListCopy[j[k].RB1.bodyID].Positions, reshape(j[k].pos1,(1,3)))
        extFListCopy[j[k].RB1.bodyID].Torques =
        vcat(extFListCopy[j[k].RB1.bodyID].Torques, reshape(ForceJoints[4:6,2*k-1],(1,3)))

        extFListCopy[j[k].RB2.bodyID].Forces =
        vcat(extFListCopy[j[k].RB2.bodyID].Forces, reshape(ForceJoints[1:3,2*k],(1,3)))
        extFListCopy[j[k].RB2.bodyID].Positions =
        vcat(extFListCopy[j[k].RB2.bodyID].Positions, reshape(j[k].pos2,(1,3)))
        extFListCopy[j[k].RB2.bodyID].Torques =
        vcat(extFListCopy[j[k].RB2.bodyID].Torques, reshape(ForceJoints[4:6,2*k],(1,3)))
    end
    return extFListCopy
end
