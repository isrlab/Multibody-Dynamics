# To implement joints between two rigid bodies.
# Location of joint specified initially.

# include("RigidBody.jl")
# include("OrientationConversion.jl")

using LinearAlgebra

mutable struct Joint
    RB1::RigidBody
    pos1::Vector{Float64} # Coordinates of joint in body-fixed frame of body 1

    RB2::RigidBody
    pos2::Vector{Float64} # Coordinates of joint in body-fixed frame of body 2

    type::String # "Revolute", "Revolute2", "Spherical", "Weld", "Free", "Spring"

    # Revolute Joint
    axis::Vector{Float64} # Axis for revolute joint

    # Spring
    k::Float64
    restLen::Float64

    # Joint Actuation: Defined in the body frame of the parent link
    jointF::Vector{Float64} # Length: 6 [F;τ]

    function Joint(body1::RigidBody, body2::RigidBody, rj1::Vector{Float64},
        rj2::Vector{Float64}; type::String="Free",
        axis::Vector{Float64}=zeros(3), k::Float64=0.0, rL::Float64=0.0,
        jointForce::Vector{Float64} = zeros(3), jointTorque::Vector{Float64} = zeros(3))
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
