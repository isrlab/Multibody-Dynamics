# To implement joints between two rigid bodies.
# Location of joint specified initially.

include("RigidBody.jl")
include("OrientationConversion.jl")

using LinearAlgebra

mutable struct Joint
    RB1::RigidBody
    pos1::Vector{Float64} # Coordinates of joint in body-fixed frame of body 1

    RB2::RigidBody
    pos2::Vector{Float64} # Coosrdinates of joint in body-fixed frame of body 2

    type::String # "Revolute", "Spherical", "Weld"

    axis::Vector{Float64}

    function Joint(body1::RigidBody, body2::RigidBody, rj1::Vector{Float64}, rj2::Vector{Float64}, type::String, axis::Vector{Float64})
        if type != "Revolute" || type != "Spherical" || type != "Weld"
            error("Unknown joint specified.")
        end

        this = new()
        this.RB1 = body1
        this.RB2 = body2
        this.pos1 = rj1
        this.pos2 = rj2
        this.type = type
        if type == "Revolute"
            this.axis = axis
        else
            this.axis = Vector{Float64}(undef,3,1)
        end
        return this
    end
end

# mutable struct RevJoint
#     RB1::RigidBody
#     pos1::Vector{Float64} # Coordinates of joint in body-fixed frame of body 1
#
#     RB2::RigidBody
#     pos2::Vector{Float64} # Coordinates of joint in body-fixed frame of body 2
#
#     axis::Vector{Float64} # Unit vector corresponding to Euler Axis of rotation
#
#     function RevJoint(body1::RigidBody, body2::RigidBody, rj1::Vector{Float64}, rj2::Vector{Float64}, axis::Vector{Float64})
#         this = new()
#         this.RB1 = body1
#         this.RB2 = body2
#         this.pos1 = rj1
#         this.pos2 = rj2
#         this.axis = axis
#         return this
#     end
# end
