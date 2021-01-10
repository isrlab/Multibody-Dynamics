# This is an implementation of rigid body dynamics
# RigidBody is a mutable structure which only contains data.
# We have functions that work on this data structure, updating the values of the member variables.
# Need to define springs, dampers, constraints.

using LinearAlgebra
# include("OrientationConversion.jl");


mutable struct extForces
    Forces::Matrix{T} where T<:Real# size n x 3 (n Forces)
    Positions::Matrix{T} where T<:Real # size n x 3 (same size as Forces)
    Torques::Matrix{T} where T<:Real # size m x 3 (m Torques)
end


mutable struct RigidBody
    m::T where T<:Real; # Mass of the rigid body.
    In::Matrix{T} where T<:Real; # Moment of inertia 3x3 matrix
    orientation::String; # "321", "quaternions". Only these two supported now.
    # This will determine how many states there are
    # 12 for euler angle, 13 for quaternions.
    bodyID::Integer # Body number in robot tree

    # Internal variables
    x::Vector{T} where T<:Real;    # State vector
    dcm::Matrix{T} where T<:Real;  # Inertial to body.
    invI::Matrix{T} where T<:Real; # inverse of moment of inertia, compute once during initialization.
    J::Matrix{T} where T<:Real # 4x4 Inertia with first elemtn a random positive number
    invJ::Matrix{T} where T<:Real # inverse of 4x4 inertia matrix
    ω::Vector{T} where T<:Real # Angular velocity in body frame
    function RigidBody(mass::T, Inertia::Matrix{T}, bodyID::Integer; orientation::String = "quaternions") where T<:Real
        # if orientation == "321"
        #     x0 = Vector{Float64}(undef,12); # 12 state vector
        # elseif lowercase(orientation) == "quaternions"
        x0 = Vector{T}(undef,14); # 14 state vector
        # else
        #     error("Unknown orientation specified")
        # end

        dcm = Matrix{T}(undef,3,3);

        if det(Inertia) != 0
            invI = inv(Inertia)
        else
            invI = Matrix{T}(I,3,3)
        end

        J = [rand(1) zeros(1,3)
             zeros(3,1) Inertia]
        invJ = [1/J[1,1] zeros(1,3)
                zeros(3,1) invI]
        xdot = function times2(in::T)::T
            return(2*in);
        end
        this = new();
        this.m = mass;
        this.In = Inertia;
        this.bodyID = bodyID
        this.orientation = lowercase(orientation);
        this.dcm = dcm
        this.x = x0;
        this.invI = invI;
        this.J = J
        this.invJ = invJ
        this.ω = Vector{T}(undef,3)
        return this;
    end
end

# mutable struct RigidBody
#     m::Real; # Mass of the rigid body.
#     In::Matrix{Real}# where T<:Real; # Moment of inertia 3x3 matrix
#     orientation::String; # "321", "quaternions". Only these two supported now.
#     # This will determine how many states there are
#     # 12 for euler angle, 13 for quaternions.
#     bodyID::Integer # Body number in robot tree
#
#     # Internal variables
#     x::Vector{Real}# where T<:Real;    # State vector
#     dcm::Matrix{Real}# where T<:Real;  # Inertial to body.
#     invI::Matrix{Real}# where T<:Real; # inverse of moment of inertia, compute once during initialization.
#     J::Matrix{Real}# where T<:Real # 4x4 Inertia with first elemtn a random positive number
#     invJ::Matrix{Real}# where T<:Real # inverse of 4x4 inertia matrix
#     ω::Vector{Real}# where T<:Real # Angular velocity in body frame
#     function RigidBody(mass::Real, Inertia::Matrix{Real}, bodyID::Integer; orientation::String = "quaternions") #where T<:Real
#         # if orientation == "321"
#         #     x0 = Vector{Float64}(undef,12); # 12 state vector
#         # elseif lowercase(orientation) == "quaternions"
#         x0 = Vector{Real}(undef,14); # 14 state vector
#         # else
#         #     error("Unknown orientation specified")
#         # end
#
#         dcm = Matrix{Real}(undef,3,3);
#
#         if det(Inertia) != 0
#             invI = inv(Inertia)
#         else
#             invI = Matrix{Real}(I,3,3)
#         end
#
#         J = [rand(1) zeros(Real,(1,3))
#              zeros(Real,(3,1)) Inertia]
#         invJ = [1/J[1,1] zeros(Real,(1,3))
#                 zeros(Real,(3,1)) invI]
#         xdot = function times2(in::Real)::Real
#             return(2*in);
#         end
#         this = new();
#         this.m = mass;
#         this.In = Inertia;
#         this.bodyID = bodyID
#         this.orientation = lowercase(orientation);
#         this.dcm = dcm
#         this.x = x0;
#         this.invI = invI;
#         this.J = J
#         this.invJ = invJ
#         this.ω = Vector{Real}(undef,3)
#         return this;
#     end
# end

function rbDynQuat(RB::RigidBody, unconstrainedF::Vector{Float64}, Fc::Vector{Float64})::Vector{Float64}

    # States are: 14 states
    # x = [x,y,z]:inertial, [u,v,w]:inertial, [β0,β1,β2,β3], [β̇0,β̇1,β̇2,β̇3]
    # unconstrainedF is a 7x1 vector, unique to each rigid body.
    # Fc is the  generalised force coming from constraints acting on the rigid body. Unique to each body. 7x1 vector.

    xdot = Vector{Float64}(undef,14)

    u = RB.x[8:10];   # Velocities in Inertial reference frame
    β = RB.x[4:7];  # Euler parameters
    βdot = RB.x[11:14] # Euler parameter derivatives

    TotalForce = unconstrainedF[1:3] + Fc[1:3]
    TotalMoment = unconstrainedF[4:7] + Fc[4:7]

    # # Udwadia's formulation using 14 States
    E = genE(β)
    xdot[1:3] = u
    xdot[8:10] = TotalForce/RB.m
    xdot[4:7] = βdot
    xdot[11:14] = 1/4*transpose(E)*RB.invJ*E*TotalMoment
    return xdot
end


function initialiseRigidBody!(RB::RigidBody,x0::Vector{T})::RigidBody where T<:Real
    # Initialise rigid body with x0
    RB.x = x0
    return RB
end

function updateRigidBody!(b::RigidBody,x::Vector{T})::RigidBody where T<:Real
    b.x = x
    b.dcm = quat2dcm(x[4:7])
    b.ω = angVel(x[4:7],x[11:14])
    return b
end

function InertialFrameAsRB()::RigidBody
    m = 0.0
    I = zeros(Real,(3,3))

    b = RigidBody(m,I,1)
    b.x = Vector{Real}([zeros(3);1.0;zeros(3);zeros(3);zeros(4)]);
    b.dcm = quat2dcm(b.x[4:7])
    b.ω = angVel(b.x[4:7],b.x[11:14])
    return b
end

function zeroExtForce()::extForces
    # External forces function to be changed later.
    Forces = zeros(1,3)
    Positions = zeros(1,3)
    Torques = zeros(1,3)
    return extF = extForces(Forces,Positions,Torques)
end

function zeroExtForceVec(l::Int64)::Vector{extForces}
    # l is the length of the vector
    v = Vector{extForces}(undef,l)
    for i=1:l
        v[i] = zeroExtForce()
    end
    return v
end

## Euler angle formulation (12 states)
# function rbDynEuler(xdot::Vector{Float64}, RB::RigidBody,
#     PositionList::Matrix{Float64}, ForceList::Matrix{Float64},
#     TorqueList::Matrix{Float64}, GravityInInertial::Vector{Float64})
#     # States are: 12 states
#     # x = [x,y,z]:inertial, [u,v,w]:body, [ϕ,θ,ψ], [ω1,ω2,ω3]:body
#     # Forces in ForceList are decribed in body reference, ForceList is nFx3.
#     # PositionList are positions where forces are acting relative to cm,
#     # decribed in body reference, PositionList is nFx3.
#     # Torques in TorqueList are in body reference. TorqueList is matrix nTx3
#
#     xpos = RB.x[1];  # x position in inertial frame
#     ypos = RB.x[2];  # y position in inertial frame
#     h    = RB.x[3];  # z position in inertial frame
#
#     u = RB.x[4:6];   # Velocities in body reference frame
#     ang = RB.x[7:9];     # Euler parameters
#     ω = RB.x[10:12]; # Angular velocity  in body reference frame
#
#     # DCM w.r.t Euler parameters. vb = C*vn, inertial to body.
#     ang2dcm(RB.dcm,ang,RB.orientation); # Have to code this up.
#
#     GravityInBody = RB.dcm*GravityInInertial;  # Convert to body reference
#
#     xdot[1:3] = RB.dcm'*u; # Velocities in inertial frame.
#     xdot[4:6] = sum(ForceList,dims=1)/RB.m + GravityInBody; # Euler's first law.
#     xdot[7:9] = # Code this
#
#     # ω dot equation -- Euler's second law
#     TotalMoment = sum(cross(PositionList[i,:],ForceList[i,:]) for i in 1:size(Fi)[1]) + sum(TorqueList,dims=1);
#     xdot[10:12] = RB.invI*(TotalMoment - cross(ω,RB.I*ω));
# end
