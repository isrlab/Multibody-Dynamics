# Simulation
# Assuming all bodies connected in a tree by joints


# include("RigidBody.jl")
# include("OrientationConversion.jl")
# include("Joint.jl")
# include("Force.jl")
# include("extF.jl")


using LinearAlgebra
using Revise
using DifferentialEquations
using StaticArrays
using BenchmarkTools

global GravityInInertial = [0.0,0.0,-9.806];

struct solSim
    # Structure to store results of Simulation for each rigid body
    r::Matrix{Float64}# where T<:Real
    v::Matrix{Float64}# where T<:Real
    β::Matrix{Float64}# where T<:Real
    βdot::Matrix{Float64}# where T<:Real
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
function simulate(tEnd::Float64,tInt::Float64,j::Vector{Joint};
                  g::Vector{Float64}=[0.0;0.0;-9.806]);# where T<:Real#,extFVec::Vector{extForces}=Vector{extForces}(undef,1))
    # Ellipsis (...) to facilitate supplying variable number of arguments to the function
    # Initial condition (all bodies connected)
    global GravityInInertial = g
    # global extFList = extFVec
    X0 = Vector{Float64}(j[1].RB1.x)
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end
    # println("X0 type = ", typeof(X0))
    # println("tEnd type = ", typeof(tEnd))
    # Declaring the ODE Problem as per DifferentialEquations convention
    tSpanODE = (0.0, tEnd);

    prob = ODEProblem(mainDynODEDiff!,X0,tSpanODE,j)
    # sol = solve(prob, saveat = tSpan, RK4(), reltol=1e-10, abstol=1e-10)
    # sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10)#, saveat = 0:tSpan:tEnd)
    sol = solve(prob, Tsit5(), saveat=tInt, reltol=1e-10, abstol=1e-10)
    # sol = solve(prob, DP5(), reltol=1e-10, abstol=1e-10)

    solFinal = Vector{solSim}(undef,length(j))
    for i=1:length(j)
        rSol = Matrix{Float64}(undef,(length(sol.t),3))
        vSol = Matrix{Float64}(undef,(length(sol.t),3))
        βSol = Matrix{Float64}(undef,(length(sol.t),4))
        βdotSol = Matrix{Float64}(undef,(length(sol.t),4))
        for t=1:length(sol.t)
            rSol[t,:] = (sol.u[t][14*i+1:14*i+3])
            vSol[t,:] = (sol.u[t][14*i+8:14*i+10])
            βSol[t,:] = (sol.u[t][14*i+4:14*i+7])
            βdotSol[t,:] = (sol.u[t][14*i+11:14*i+14])
        end
    solFinal[i] = solSim(rSol,vSol,βSol,βdotSol)
    end
    return (sol.t, solFinal)
end

function mainDynODEDiff!(dX,X, j::Vector{Joint}, t::Float64)# where T<:Real
    dX[:] = mainDynODEDiff(X,j,t);
end

function mainDynODEDiff(X::Vector{T}, j::Vector{Joint}, t::Float64)::Vector{T} where T<:Real
    # ODE function to be used as per DifferentialEquations convention
    # Create extForcesList storing extForces for each rigid body
    # Create ForceConstraints Array storing constraint forces acting on each rigid body
    global GravityInInertial
    # Update RigidBodies
    updateRigidBody!(j[1].RB1,X[1:14]) # update inertial frame (not moving)
    for k=1:length(j)
        updateRigidBody!(j[k].RB2,X[14*k+1:14*(k+1)])
    end

    extFListCopy = externalFTotal(t,j); # Generates all the external forces explicitly specified, through joint actions or directly.

    U = genU(extFListCopy) # Construct an input matrix from external forces and torques supplied.

    ## Inplace
    dX = zeros(length(X));
    fxdot!(dX, X, U,j,GravityInInertial);

    ## Not inplace
    # dX = fxdot(X,U,j,GravityInInertial);

    # println("t = $t")
    return dX
end

function fxdot!(dX, X::Vector{T}, U::Matrix{S}, j::Vector{Joint}, GravityInInertial::Vector{Float64}) where {T<:Real,S<:Real}
    dX[:] = fxdot(X,U,j,GravityInInertial);
end

function fxdot(X::Vector{T},U::Matrix{S},j::Vector{Joint},GravityInInertial::Vector{Float64}) where {T<:Real,S<:Real}
    nJ = length(j); # Number of joints

    unconstrF, constrF = Constraint(X, U, j, GravityInInertial)

    xdot = Vector{Union{T,S}}(undef,14*(nJ+1));
    xdot[1:14] = zeros(14);

    for k=1:nJ
        xb = X[14*k+1:14*(k+1)]

        unconstrainedF_rb = unconstrF[:,j[k].RB2.bodyID]
        Fc_rb = constrF[:,j[k].RB2.bodyID]
        TotalForce = unconstrainedF_rb[1:3] + Fc_rb[1:3]
        TotalMoment = unconstrainedF_rb[4:7] + Fc_rb[4:7]
        genF = [TotalForce; TotalMoment]; # generalized force acting on RB
        genM = Array(genMatM(X, j[k].RB2)) # generalized mass matrix of RB
        xdot[14*k+1:14*k+7] = X[14*k+8:14*k+14]#xb[8:10]
        xdot[14*k+8:14*k+14] = inv(genM)*genF;
    end

    return xdot
end

function externalFTotal(t::T, j::Vector{Joint})::Vector{extForces} where T<:Real
    # Returns externally applied forces in total, not including gravity. Includes Forces from joints and other explicitly specified external forces.
    extFList = extF(t,j)

    # Generate forces from actuated joints on each body
    # ForceJoints = Matrix{Float64}(undef,6,2*length(j))
    # for k=1:length(j)
    #     ForceJoints[:,2*k-1], ForceJoints[:,2*k] = genJointF(t,j[k])
    # end

    # Add ForceJoints to extFList
    extFListCopy = deepcopy(extFList)
    # for k=1:length(j)
    #     extFListCopy[j[k].RB1.bodyID].Forces =
    #     vcat(extFListCopy[j[k].RB1.bodyID].Forces, reshape(ForceJoints[1:3,2*k-1],(1,3)))
    #     extFListCopy[j[k].RB1.bodyID].Positions =
    #     vcat(extFListCopy[j[k].RB1.bodyID].Positions, reshape(j[k].pos1,(1,3)))
    #     extFListCopy[j[k].RB1.bodyID].Torques =
    #     vcat(extFListCopy[j[k].RB1.bodyID].Torques, reshape(ForceJoints[4:6,2*k-1],(1,3)))
    #
    #     extFListCopy[j[k].RB2.bodyID].Forces =
    #     vcat(extFListCopy[j[k].RB2.bodyID].Forces, reshape(ForceJoints[1:3,2*k],(1,3)))
    #     extFListCopy[j[k].RB2.bodyID].Positions =
    #     vcat(extFListCopy[j[k].RB2.bodyID].Positions, reshape(j[k].pos2,(1,3)))
    #     extFListCopy[j[k].RB2.bodyID].Torques =
    #     vcat(extFListCopy[j[k].RB2.bodyID].Torques, reshape(ForceJoints[4:6,2*k],(1,3)))
    # end
    return extFListCopy
end

function genU(extFList::Vector{extForces})::Matrix{Float64}
    nB = length(extFList); # number of bodies
    u = zeros(6,nB); # Each column in u contains a vector of [F;τ] acting on the corresponding rigid body
    for b=1:nB
        extF_b = extFList[b]
        u[1:3,b] = sum(extF_b.Forces,dims = 1) # Sum up forces acting on body
        for n in 1:size(extF_b.Forces,1) # Sum up moments acting on body
            u[4:6,b] += cross(extF_b.Positions[n,:], extF_b.Forces[n,:])
        end
        u[4:6,b] += sum(extF_b.Torques, dims = 1)[:] # Sum up torques acting on body
    end
    return u
end


function getXU_0(j::Vector{Joint})
    X0 = Vector{Float64}(j[1].RB1.x)
    for k = 1:length(j)
        append!(X0,j[k].RB2.x)
    end

    extFListCopy = externalFTotal(0.0,j); # Generates all the external forces explicitly specified, through joint actions or directly.

    U0 = genU(extFListCopy) # Construct an input matrix from external forces and torques supplied.

    return X0, U0
end

function checkRevJoint(solQuad::solSim, solProp1::solSim, rjCube1::Array{Float64}, rjProp1::Array{Float64})
    ## Function that outputs relevant quantities constrained under a rotational joint
    tLen = size(solQuad.r,1)
    ωCube = Matrix{Float64}(undef,tLen,3);
    ωProp1 = Matrix{Float64}(undef,tLen,3);
    ωProp1InCube = Matrix{Float64}(undef,tLen,3);
    ωRel = Matrix{Float64}(undef,tLen,3);
    jointLoc1 = Matrix{Float64}(undef,length(tSim),3);
    for i=1:tLen
        dcm = quat2dcm(solQuad.β[i,:])
        dcm2 = quat2dcm(solProp1.β[i,:])
        ωCube[i,:] = angVel(solQuad.β[i,:],solQuad.βdot[i,:])
        ωProp1[i,:] = angVel(solProp1.β[i,:],solProp1.βdot[i,:])
        ωProp1InCube[i,:] = quat2dcm(solQuad.β[i,:])*transpose(quat2dcm(solProp1.β[i,:]))*ωProp1[i,:]
        ωRel[i,:] = ωProp1InCube[i,:] - ωCube[i,:]
        jointLoc1[i,:] = solQuad.r[i,:] + transpose(dcm)*rjCube1 - solProp1.r[i,:] - transpose(dcm2)*rjProp1;
        # jointLoc1[i,:] = solProp1.r[i,:] + transpose(dcm2)*rjProp1;
    end
    return ωCube, ωProp1, ωRel, jointLoc1
end

function checkRevJointIn(solPend1::solSim, rj2::Array{Float64})
    ## Function that outputs relevant quantities constrained under a rotational joint when first body is inertial frame
    tLen = size(solPend1.r,1)
    jointLoc1 = Matrix{Float64}(undef,tLen,3)
    ωSol = Matrix{Float64}(undef,tLen,3)
    for i=1:tLen
        tmp_dcm = quat2dcm(solPend1.β[i,:])
        ωSol[i,:] = angVel(solPend1.β[i,:],solPend1.βdot[i,:])
        jointLoc1[i,:] = solPend1.r[i,:] + transpose(tmp_dcm)*rj2
    end
    return ωSol, jointLoc1
end

function joinBodies!(j1::Vector{Joint}, j2::Vector{Joint}, b1_id::Integer, b2_id::Integer, rj1::Vector{T}, rj2::Vector{T}; type::String="Weld", axis::Vector{T}=[0.0;0.0;1.0], k::T=0.0, rL::T=0.0, jointForce::Vector{T} = zeros(T,3), jointTorque::Vector{T} = zeros(T,3))  where T<:Real
    # join the bodies b1 given by b1_id in j1 and b2_id in j2 at rj1 and rj2 respectively
    # Create new joint as a weld joint, unless specified
    len_j1 = length(j1); # old length of j1
    RB1 = body_in_jointVec(j1,b1_id);
    RB2 = body_in_jointVec(j2,b2_id);

    newJ = Joint(RB1, RB2, rj1, rj2, type = type, axis = axis, k = k, rL = rL, jointForce = jointForce, jointTorque = jointTorque);
    push!(j1,newJ)

     # Update bodyIDs of all bodies present in joint tree j2
    for i=1:length(j2)
        j2[i].RB2.bodyID += len_j1;
        push!(j1,j2[i])
    end
    filterInertialJoint!(j1); # This removes the connection of the body in joint tree j2 that was connected to the inertial frame (removes loop inside joint trees)
end
