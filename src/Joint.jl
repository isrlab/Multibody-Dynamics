# To implement joints between two rigid bodies.
# Location of joint specified initially.

# include("RigidBody.jl")
# include("OrientationConversion.jl")

using LinearAlgebra

mutable struct Joint
    RB1::RigidBody
    pos1::Vector{T} where T<:Real # Coordinates of joint in body-fixed frame of body 1

    RB2::RigidBody
    pos2::Vector{T} where T<:Real# Coordinates of joint in body-fixed frame of body 2

    type::String # "Revolute", "Revolute2", "Spherical", "Weld", "Free", "Spring"

    # Revolute Joint
    axis::Vector{T} where T<:Real # Axis for revolute joint

    # Spring
    k::T where T<:Real
    restLen::T where T<:Real

    # Joint Actuation: Defined in the body frame of the parent link
    jointF::Vector{T} where T<:Real # Length: 6 [F;τ]

    function Joint(body1::RigidBody, body2::RigidBody, rj1::Vector{T},
        rj2::Vector{T}; type::String="Free",
        axis::Vector{T}=[0.0;0.0;1.0], k::T=0.0, rL::T=0.0,
        jointForce::Vector{T} = zeros(T,3), jointTorque::Vector{T} = zeros(T,3)) where T<:Real
        # Default values for type, axis, k, and restLen provided.
        allowedTypes = ["Revolute", "Revolute2", "Spherical", "Weld", "Free", "Spring"]
        if !in(type,allowedTypes)
            error("Unknown joint specified.")
        end

        this = new()
        this.RB1 = body1
        this.RB2 = body2
        this.pos1 = rj1
        this.pos2 = rj2
        this.type = type
        this.axis = axis
        this.k = k
        this.restLen = rL
        this.jointF = [jointForce;jointTorque]

        return this
    end
end

# Number of constraint equations for all kinds of joints
JointDict = Dict(["Revolute"=>7, "Revolute2"=>6, "Spherical" => 5, "Weld"=>9, "Spring" => 2, "Free" => 2]);
JointDictIn = Dict(["Revolute"=>6, "Revolute2"=>5, "Spherical" => 4, "Weld"=>8 , "Spring" => 1, "Free" => 1]); # When the first body is the inertial frame

function body_in_jointVec(jVec::Vector{Joint}, b_id::Int64)
# returns the joint id which contains b_id in the vector of joints
    # j = 0;
    # j_i = 0;
    RB = InertialFrameAsRB();
    for i=1:length(jVec)
        if (b_id == jVec[i].RB1.bodyID)
            # j_i = i;
            RB = jVec[i].RB1;
            break;
        elseif (b_id ==  jVec[i].RB2.bodyID)
            # j_i = i;
            # j = 2;
            RB = jVec[i].RB2;
            break;
        else continue;
        end
    end
    return RB;
    # j_i is the position of the first joint in the vector that contains corresponding body
    # j is either 1 or 2, i.e., which body of the corresponding joint
end

function filterInertialJoint!(jVec::Vector{Joint})
    # Remove the inertial body/frame occurring twice in a joint tree
    for i=2:length(jVec)
        if (jVec[i].RB1.m == 0 || jVec[i].RB2.m == 0)
            deleteat!(jVec,i);
            break;
        end
        # RB1 = jVec[i].RB1; RB2 = jVec[i].RB2;
        # for j=1:length(jVec)
        #
        # end
    end
end
