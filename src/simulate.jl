# Simulation
# Tentatively to simulate two rigid bodies constrained by one revolute joint


include("RigidBody.jl")
include("OrientationConversion.jl")
include("Joint.jl")
include("Force.jl")
include("extF.jl")


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

"""
simulate(tEnd,tSpan,j;g,extFVec)

tEnd: Simulation End time. (Float64) \\
tSpan: Result saved at every specified tSpan. (Float64) \\
j: Tuple of Joints. Length need not be specified. \\
g: Acceleration due to gravity. Default value: [0.0;0.0;-9.806] (Vector{Float64}) \\
extFVec: Vector of extForces acting on each rigid body. (Vector{extForces}) \\

Tolerance of ode solver set. RelTol = 1e-10, AbsTol = 1e-10. \\
Tsit5() solver used. \\
Returns a tuple containing tSim and vector of solutions in solSim form.
"""
function simulate(tEnd::Float64,tSpan::Float64,j::Joint...;
                  g::MArray{Tuple{3},Float64,1,3}=[0.0;0.0;-9.806])#,extFVec::Vector{extForces}=Vector{extForces}(undef,1))
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function
    # Initial condition (all bodies connected)
    global GravityInInertial = g
    # global extFList = extFVec
    X0 = j[1].RB1.x
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end
    # Declaring the ODE Problem as per DifferentialEquations convention
    prob = ODEProblem(mainDynODE!,X0,(0.0,tEnd),j)
    sol = solve(prob,Tsit5(),reltol=1e-10,abstol=1e-10)
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


"""
mainDynODE!(dX,X,j,t) \\
X:: State vector of system at time t. \\
j:: Tuple of Joints. (Length not specified.) \\

Generates the constraint forces acting on each rigid body present in the system at time t. \\
dX = f(X,t)
"""
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
    @show t

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
