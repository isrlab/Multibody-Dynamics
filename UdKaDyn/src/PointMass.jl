# This is an implementation of point mass dynamics
# PointMass is a mutable structure which only contains data.
# We have functions that work on this data structure, updating the values of the member variables.

using LinearAlgebra

mutable struct PointMass
    m::Float64; # Mass of the rigid body.

    # Internal variables
    x::Vector{Float64};    # State vector

    function PointMass(mass::Float64)
        this = new();
        this.m = mass;
        this.In = Inertia;
        return this;
    end
end

function pmDyn(PM::PointMass,
    extF::extForces, Fc::Vector{Float64}, GravityInInertial::Vector{Float64})

    # States are: 6
    # x = [x,y,z]:inertial, [u,v,w]:inertial
    # extF and Fc are described in the inertial frame.

    xdot = Vector{Float64}(undef,6)

    unconstrainedF = genExtF_PM(PM,extF,GravityInInertial)
    TotalForce = unconstrainedF[1:3] + Fc[1:3]

    # Udwadia's formulation using 6 States
    xdot[1:3] = x[4:6]
    xdot[4:6] = TotalForce/PM.m
    return xdot
end


function initialisePointMass(p::PointMass,x0::Vector{Float64})
    # Initialise point mass with x0
    p.x = x0
    return RB
end

function updatePointMass(p::PointMass,x::Vector{Float64})
    p.x = x
    return p
end
