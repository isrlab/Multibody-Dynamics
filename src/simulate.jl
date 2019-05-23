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

global GravityInInertial = Vector{Float64}(undef,3)
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
                  g::Vector{Float64}=[0.0;0.0;-9.806])#,extFVec::Vector{extForces}=Vector{extForces}(undef,1))
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function
    # Initial condition (all bodies connected)
    global GravityInInertial = g
    # global extFList = extFVec
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


"""
mainDynODE(X,j,t) \\
X:: State vector of system at time t. \\
j:: Tuple of Joints. (Length not specified.) \\

Generates the constraint forces acting on each rigid body present in the system at time t. \\
Supplies these forces to mainDyn function.
"""
function mainDynODE(X::Vector{Float64},j::Tuple{Vararg{Joint}},t::Float64)
    # ODE function to be used as per DifferentialEquations covnention
    # Create extForcesList storing extForces for each rigid body
    # Create ForceConstraints Array storing constraint forces acting on each rigid body
    global GravityInInertial
    extFList = extF(t,j...)
    # Update RigidBodies
    updateRigidBody!(j[1].RB1,X[1:14])
    for k=1:length(j)
        updateRigidBody!(j[k].RB2,X[14*k+1:14*(k+1)])
    end
    @show t

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

    unconstrF, ForceConstr = Constraint(j,extFListCopy,GravityInInertial)

    ## Generate constraint Forces for each body using UK Formulation
    # ForceConstr = zeros(7,length(j)+1)
    # ForceConstr[:,2] = ForceCon(j[1],extFListCopy[1],extFListCopy[2],
    #                    ForceConstr[:,1],ForceConstr[:,2],GravityInInertial)
    # if length(j) > 1
    #     for k = 2:length(j)
    #         FcAug = ForceCon(j[k],extFListCopy[j[k].RB1.bodyID],
    #         extFListCopy[j[k].RB2.bodyID],ForceConstr[:,j[k].RB1.bodyID],
    #         ForceConstr[:,j[k].RB2.bodyID],GravityInInertial)
    #
    #         ForceConstr[:,j[k].RB1.bodyID] =
    #         ForceConstr[:,j[k].RB1.bodyID] + FcAug[1:7]
    #         # ForceCon(j[k],extFListCopy[j[k].RB1.bodyID],
    #         # extFListCopy[j[k].RB2.bodyID],ForceConstr[:,j[k].RB1.bodyID],
    #         # ForceConstr[:,j[k].RB2.bodyID],GravityInInertial)[1:7]
    #
    #         ForceConstr[:,j[k].RB2.bodyID] =
    #         ForceConstr[:,j[k].RB2.bodyID] + FcAug[8:14]
    #         # ForceCon(j[k],extFListCopy[j[k].RB1.bodyID],
    #         # extFListCopy[j[k].RB2.bodyID],ForceConstr[:,j[k].RB1.bodyID],
    #         # ForceConstr[:,j[k].RB2.bodyID],GravityInInertial)[8:14]
    #
    #         # if k==3
    #         #     # a = ForceCon(j[k],extFListCopy[j[k].RB1.bodyID],extFListCopy[j[k].RB2.bodyID],GravityInInertial)[1:7]
    #         #     @show a
    #         #     # @show ForceConstr[:,3]
    #         # end
    #         # if k==2
    #         #     # @show extFListCopy[3]
    #         #     @show ForceConstr[:,3]
    #         # end
    #     end
    # end

    dX = mainDyn(X,j,unconstrF,ForceConstr, GravityInInertial)
    return dX
end

"""
mainDyn(Q,j,extFList,ForceConstr, GravityInInertial) \\

Q:: State Vector at time t.\\
j:: Tuple of Joints \\
extFList:: Vector of extForces acting on each rigid body \\
ForceConstr:: Constraint Forces generated using U-K formulation, acting on each rigid body.

Main function for solving the ODE. \\
Output: dQ = f(Q,t).
"""
function mainDyn(Q::Vector{Float64},j::Tuple{Vararg{Joint}},
    unconstrF::Matrix{Float64}, ForceConstr::Matrix{Float64},
    GravityInInertial::Vector{Float64})

    dQ = Vector{Float64}(undef,(length(j)+1)*14)

    # First body always the inertial frame
    dQ[1:14] = zeros(14)
    for k=1:length(j)
        # rbDynQuat!(dQ[14*k+1:14*(k+1)],j[k].RB2,extFList[k+1],ForceConstr[:,k+1],GravityInInertial)
        dQ[14*k+1:14*(k+1)] = rbDynQuat(j[k].RB2,unconstrF[:,j[k].RB2.bodyID],
                              ForceConstr[:,j[k].RB2.bodyID],GravityInInertial)

    end
    return dQ
end

# function getExternalForce(b::RigidBody)
#     # External forces function to be changed later.
#     Forces = zeros(1,3)
#     Positions = zeros(1,3)
#     Torques = zeros(1,3)
#     return extF = extForces(Forces,Positions,Torques)
# end
