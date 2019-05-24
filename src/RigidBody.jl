# This is an implementation of rigid body dynamics
# RigidBody is a mutable structure which only contains data.
# We have functions that work on this data structure, updating the values of the member variables.
# Need to define springs, dampers, constraints.

using LinearAlgebra
# include("OrientationConversion.jl");


mutable struct extForces
    Forces::Matrix{Float64} # size n x 3 (n Forces)
    Positions::Matrix{Float64} # size n x 3 (same size as Forces)
    Torques::Matrix{Float64} # size m x 3 (m Torques)
end


mutable struct RigidBody
    m::Float64; # Mass of the rigid body.
    In::Matrix{Float64}; # Moment of inertia 3x3 matrix
    orientation::String; # "321", "quaternions". Only these two supported now.
    # This will determine how many states there are
    # 12 for euler angle, 13 for quartenions.
    bodyID::Integer # Body number in robot tree

    # Internal variables
    x::Vector{Float64};    # State vector
    dcm::Matrix{Float64};  # Inertial to body.
    invI::Matrix{Float64}; # inverse of moment of inertia, compute once during initialization.
    J::Matrix{Float64} # 4x4 Inertia with first elemtn a random positive number
    invJ::Matrix{Float64} # inverse of 4x4 inertia matrix
    ω::Vector{Float64} # Angular velocity in body frame
    function RigidBody(mass::Float64, Inertia::Matrix{Float64}, bodyID::Integer; orientation::String = "quaternions")
        # if orientation == "321"
        #     x0 = Vector{Float64}(undef,12); # 12 state vector
        # elseif lowercase(orientation) == "quaternions"
        x0 = Vector{Float64}(undef,14); # 14 state vector
        # else
        #     error("Unknown orientation specified")
        # end

        dcm = Matrix{Float64}(undef,3,3);

        if det(Inertia) != 0
            invI = inv(Inertia)
        else
            invI = Matrix{Float64}(I,3,3)
        end

        J = [rand(1) zeros(1,3)
             zeros(3,1) Inertia]
        invJ = [1/J[1,1] zeros(1,3)
                zeros(3,1) invI]
        xdot = function times2(in::Float64)::Float64
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
        this.ω = Vector{Float64}(undef,3)
        return this;
    end
end

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

function rbDynQuat(RB::RigidBody, unconstrainedF::Vector{Float64}, Fc::Vector{Float64})

    # States are: 14 states
    # x = [x,y,z]:inertial, [u,v,w]:inertial, [β0,β1,β2,β3], [β̇0,β̇1,β̇2,β̇3]
    # unconstrainedF is a 7x1 vector, unique to each rigid body.
    # Fc is the force generated from constraints acting on the rigid body. Unique to each body. 7x1 vector.

    xdot = Vector{Float64}(undef,14)

    xpos = RB.x[1];  # x position in inertial frame
    ypos = RB.x[2];  # y position in inertial frame
    h    = RB.x[3];  # z position in inertial frame

    u = RB.x[8:10];   # Velocities in Inertial reference frame
    β = RB.x[4:7];  # Euler parameters
    βdot = RB.x[11:14] # Euler parameter derivatives
    RB.ω = angVel(β,βdot)
    β0 = β[1];
    β1 = β[2];
    β2 = β[3];
    β3 = β[4];

    # DCM w.r.t Euler parameters. vb = C*vn, inertial to body.
    RB.dcm = quat2dcm(β);

    # unconstrainedF = genExtF(RB,extF,GravityInInertial)
    TotalForce = unconstrainedF[1:3] + Fc[1:3]
    TotalMoment = unconstrainedF[4:7] + Fc[4:7]

    # Udwadia's formulation using 14 States
    E = genE(β)
    xdot[1:3] = u
    xdot[8:10] = TotalForce/RB.m
    xdot[4:7] = βdot
    xdot[11:14] = 1/4*transpose(E)*RB.invJ*E*TotalMoment
    return xdot
end


function initialiseRigidBody!(RB::RigidBody,x0::Vector{Float64})
    # Initialise rigid body with x0
    RB.x = x0
    return RB
end

function updateRigidBody!(b::RigidBody,x::Vector{Float64})
    b.x = x
    b.dcm = quat2dcm(x[4:7])
    b.ω = angVel(x[4:7],x[11:14])
    return b
end

function InertialFrameAsRB()
    m = 0.0
    I = zeros(3,3)

    b = RigidBody(m,I,1)
    b.x = [zeros(3);1;zeros(3);zeros(3);zeros(4)]
    b.dcm = quat2dcm(b.x[4:7])
    b.ω = angVel(b.x[4:7],b.x[11:14])
    return b
end

function zeroExtForce()
    # External forces function to be changed later.
    Forces = zeros(1,3)
    Positions = zeros(1,3)
    Torques = zeros(1,3)
    return extF = extForces(Forces,Positions,Torques)
end

function zeroExtForceVec(l::Int64)
    # l is the length of the vector
    v = Vector{extForces}(undef,l)
    for i=1:l
        v[i] = zeroExtForce()
    end
    return v
end
