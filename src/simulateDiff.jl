# Simulation
# Assuming all bodies connected in a tree by joints


include("RigidBody.jl")
include("OrientationConversion.jl")
include("Joint.jl")
include("ForceDiff.jl")
include("extF.jl")


using LinearAlgebra
using Revise
using DifferentialEquations
using StaticArrays

global GravityInInertial = MVector{3,Real}([0.0,0.0,-9.806])

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
function simulate(tEnd::Float64,tInt::T,j::Tuple{Vararg{Joint}};
                  g::MArray{Tuple{3},Real,1,3}=[0.0;0.0;-9.806]) where T<:Real#,extFVec::Vector{extForces}=Vector{extForces}(undef,1))
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
    tSpanODE = (0.0, tEnd)
    # println("tSpanODE type = ", typeof(tSpanODE))
    prob = ODEProblem(mainDynODE,X0,tSpanODE,j)
    # sol = solve(prob, saveat = tSpan, RK4(), reltol=1e-10, abstol=1e-10)
    # sol = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10)#, saveat = 0:tSpan:tEnd)
    sol = solve(prob, DP5(), reltol=1e-10, abstol=1e-10)

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
function mainDynODE(X::Vector{T}, j::Tuple{Vararg{Joint}}, t::Float64)::Vector{T} where T<:Real

    # println("\ntypeof(X) = ", typeof(X))
    # nB = j[end].RB2.bodyID - 1; # Number of bodies, discounting inertial frame
    #
    # nC = 0 # Number of constraint equations formed
    # rows = Int64[]
    # for i=1:nJ
    #     nC += get(JointDict, j[i].type, 1)
    #     push!(rows, nC)
    # end
    # X = deepcopy(X0);t = 0.0;
    # X = Vector{T}(undef,length(nJ+1));t = T(0.0);

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

    # constrF = Constraint(X, j, extFListCopy, GravityInInertial)
    # constrFJac = ForwardDiff.jacobian(x -> Constraint(x,j,extFListCopy, GravityInInertial), X) # works if Constraint fn returns a vector
##
    # function rbDyn(xb::Vector{T},unconstrainedF,Fc, rb)::Vector{Float64} where T<:Real
    #     f = zeros(14,1) #where T<:Real
    #     TotalForce = unconstrainedF[1:3] + Fc[1:3]
    #     TotalMoment = unconstrainedF[4:7] + Fc[4:7]
    #     f[1:3] = xb[8:10]
    #     f[8:10] = TotalForce/rb.m
    #     f[4:7] = xb[11:14]
    #     E = genE(xb[4:7])
    #     f[11:14] = 1/4*transpose(E)*rb.invJ*E*TotalMoment
    #     return f
    # end
    # clearconsole()


    dX = fxdot(X,j,extFListCopy,GravityInInertial);#[3];
    println("t = $t")
    # println("Here")
    # sleep(1000)
    # println("uF = ", uF)
    # println("cF = ", cF)
    # println(" X = ", X[29:end])
    # println("dX = ", dX[15:end])
    #
    # jb = ForwardDiff.jacobian(z -> fxdot([X[1:14];z],j,extFListCopy,GravityInInertial)[3][15:end],X[15:end])
    # jb = ForwardDiff.jacobian(z -> fxdot(z,j,extFListCopy,GravityInInertial),X)
    # println("jb = ", jb[15:end, 15:end])
    # sleep(1000)
##
    # unconstrF is the matrix of the unconstrained Forces acting on the robot.
    # constrF is the matrix of the constraint Forces acting on the system developed using UK formulation

    # dX[1:14] = zeros(14) # For the inertial frame (no motion)
    # for k=1:length(j)
    #     dX[14*k+1:14*(k+1)] = rbDynQuat(j[k].RB2, unconstrF[:,j[k].RB2.bodyID], constrF[:,j[k].RB2.bodyID])
    # end
    return dX
end

function externalFTotal(t::T, j::Tuple{Vararg{Joint}})::Vector{extForces} where T<:Real
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

function fxdot(X::Vector{T},j,extFListCopy,GravityInInertial) where T<:Real
    # println("typeof(X) = ", typeof(X));
    nJ = length(j); # Number of joints
    unconstrF, constrF = Constraint(X, j, extFListCopy, GravityInInertial)
    # println("unconstrF[:,2]= ", unconstrF[:,2])
    # println("constrF[:,2]= ", constrF[:,2])
    # sleep(1000);
    # jb = ForwardDiff.jacobian(z -> Constraint([X[1:14];z], j, extFListCopy, GravityInInertial)[2],X[15:end])
    # println("jb = ", jb)
    # sleep(1000);
    xdot = zeros(T,(14*(nJ+1)))# where T<:Real
    # xdot = zeros(3)
    xdot[1:14] = zeros(T,14)

    # println(ForwardDiff.jacobian(x -> Constraint(x,j,extFListCopy,GravityInInertial)[2][:],X))
    for k=1:nJ
        xb = X[14*k+1:14*(k+1)]

        unconstrainedF_rb = unconstrF[:,j[k].RB2.bodyID]
        Fc_rb = constrF[:,j[k].RB2.bodyID]
        TotalForce = unconstrainedF_rb[1:3] + Fc_rb[1:3]
        TotalMoment = unconstrainedF_rb[4:7] + Fc_rb[4:7]

        # xdot = xb[8:10]
        xdot[14*k+1:14*k+3] = X[14*k+8:14*k+10]#xb[8:10]
        xdot[14*k+4:14*k+7] = X[14*k+11:14*k+14]
        xdot[14*k+8:14*k+10] = TotalForce/j[k].RB2.m
        E = genE(X[14*k+4:14*k+7])
        xdot[14*k+11:14*k+14] = 1/4*transpose(E)*j[k].RB2.invJ*E*TotalMoment
        # xdot[14*k+1:14*(k+1)] = rbDyn(xb,unconstrF[:,j[k].RB2.bodyID], constrF[:,j[k].RB2.bodyID],j[k].RB2)
    end

    # return unconstrF, constrF, xdot
    return xdot
end

function checkRevJoint(solQuad::solSim, solProp1::solSim, rjCube1::Array{Float64}, rjProp1::Array{Float64})
    ## Function that outputs relevant quantities constrained under a revolute joint
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
