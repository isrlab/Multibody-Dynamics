# This is an implementation of rigid body dynamics
# RigidBody is a mutable structure which only contains data.
# We have functions that work on this data structure, updating the values of the member variables.
# Need to define springs, dampers, constraints.

using LinearAlgebra
include("OrientationConversion.jl");

mutable struct extForces
    Forces::Matrix{Float64}
    Positions::Matrix{Float64}
    Torques::Matrix{Float64}
end


mutable struct RigidBody
    m::Float64; # Mass of the rigid body.
    I::Matrix{Float64}; # Moment of inertia 3x3 matrix
    orientation::String; # "321", "quaternions". Only these two supported now.
    # This will determine how many states there are
    # 12 for euler angle, 13 for quartenions.

    # Internal variables
    x::Vector{Float64};    # State vector
    dcm::Matrix{Float64};  # Inertial to body.
    invI::Matrix{Float64}; # inverse of moment of inertia, compute once during initialization.
    J::Matrix{Float64} # 4x4 Inertia with first elemtn a random positive number
    ω::Vector{Float64} # Angular velocity in body frame
    function RigidBody(mass::Float64, Inertia::Matrix{Float64}, orientation::String)
        if orientation == "321"
            x0 = Vector{Float64}(undef,12); # 12 state vector
        elseif lowercase(orientation) == "quaternions"
            x0 = Vector{Float64}(undef,14); # 14 state vector
        else
            error("Unknown orientation specified")
        end

        dcm = Matrix{Float64}(undef,3,3);
        invI = inv(Inertia);
        J = [rand(1) zeros(1,3)
             zeros(3,1) Inertia]
        xdot = function times2(in::Float64)::Float64
            return(2*in);
        end
        this = new();
        this.m = mass;
        this.I = Inertia;
        this.orientation = lowercase(orientation);
        this.dcm = dcm
        this.x = x0;
        this.invI = inv(Inertia);
        this.J = J
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

function rbDynQuat!(xdot::Vector{Float64}, RB::RigidBody,
    extF::extForces, Fc::Vector{Float64}, GravityInInertial::Vector{Float64})

    # States are: 14 states
    # x = [x,y,z]:inertial, [u,v,w]:inertial, [β0,β1,β2,β3], [β̇0,β̇1,β̇2,β̇3]
    # Forces in ForceList are decribed in body reference, ForceList is nFx3.
    # PositionList are positions where forces are acting relative to cm,
    # decribed in body reference, PositionList is nFx3.
    # Torques in TorqueList are in body reference. TorqueList is matrix nTx3
    # Inertia is a data structure: Inertia.I and Inertia.invI.

    xpos = RB.x[1];  # x position in inertial frame
    ypos = RB.x[2];  # y position in inertial frame
    h    = RB.x[3];  # z position in inertial frame

    u = RB.x[8:10];   # Velocities in Inertial reference frame
    β = RB.x[4:7];  # Euler parameters
    βdot = RB.x[11:14] # Euler parameter derivatives
    angVel(RB.ω,β,βdot)
    β0 = β[1];
    β1 = β[2];
    β2 = β[3];
    β3 = β[4];

    ForceList = extF.Forces
    PositionList = extF.Positions
    TorqueList = extF.Torques

    # DCM w.r.t Euler parameters. vb = C*vn, inertial to body.
    quat2dcm(RB.dcm,β);

    GravityInBody = RB.dcm*GravityInInertial; # Convert to body reference

    # xdot[1:3] = RB.dcm'*u; # Velocities in inertial frame.
    # xdot[4:6] = sum(ForceList,dims=1)/RB.m + GravityInBody; # Euler's first law.
    # betadot equations
    # Ω = [0    -ω[1] -ω[2]  -ω[3];
    #      ω[1]  0     ω[2]   ω[3];
    #      ω[2] -ω[3]  0      ω[1];
    #      ω[3]  ω[2] -ω[1]    0;
    #     ];
    # xdot[7:10] = Ω*β; # Not accurate -- need Udwadia's formulation.
    # ω dot equation -- Euler's second law
    TotalForce = sum(ForceList,dims=1) + m*GravityInInertial + Fc[1:3]
    TotalMoment = sum(cross(PositionList[i,:],ForceList[i,:]) for i in 1:size(ForcesList)[1]) + sum(TorqueList,dims=1) + Fc[4:end]
    # xdot[11:13] = RB.invI*(TotalMoment - cross(ω,RB.I*ω));

    # Udwadia's formulation using 14 States
    E1 = [-β1  β0  β3 -β2
          -β2 -β3  β0  β1
          -β3  β2 -β1  β0]
    xdot[1:3] = u
    xdot[8:10] = transpose(RB.dcm)*sum(ForceList,dims=1)/RB.m + GravityInInertial
    xdot[4:7] = βdot
    xdot[11:14] = -0.5*transpose(E1)*RB.invI*skewX(RB.ω)*RB.I*RB.ω - (transpose(βdot)*βdot)*β + transpose(E1)*RB.invI*E1*TotalMoment/4
end

function initialiseRigidBody(RB::RigidBody,x0::Vector{Float64})
    # Initialise rigid body with x0
    RB.x = x0
    return RB
end
