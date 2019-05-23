#
include("RigidBody.jl")
include("Joint.jl")
include("OrientationConversion.jl")

function extF(t::Float64,j::Joint...)
    # Function to generate external forces
    # that may or may not be functions of time
    extFList = zeroExtForceVec(length(j)+1)

    # First body always the inertial frame
    # extFList[1] = zeroExtForce()

    # Pendulum Test
    # extFList[2] = zeroExtForce()

    # Joint Test
    # extFList[2] = extForces(zeros(1,3), zeros(1,3), [0.0 0.0 0.01])

    # Cube Prop Test
    # extFList[2] = extForces(zeros(1,3), zeros(1,3), [0.0 0.0 0.001])

    # CoAxCop
    # extFList[2] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.0])
    # extFList[3] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.0])
    # extFList[4] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.0])
    # extFList[5] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.0])

    # Gimbal test
    # extFList[2] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.0])
    # extFList[3] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.0])

    # ModRob
    # extFList[4] = extForces(zeros(1,3),zeros(1,3),[0.0 0.0 0.001])
    return extFList
end

function genJointF(t::Float64,j::Joint)
    # For actuated joints.
    # j.jointF gives us the actuation force and torque applied ...
    # ... on the child link in the body frame of the parent link
    β1 = j.RB1.x[4:7]
    β2 = j.RB2.x[4:7]

    Fb1 = -j.jointF # Joint Force acting on parent link
    Fb2 = Vector{Float64}(undef,6) # Force acting on child link
    Fb2[1:3] = quat2dcm(β2)*transpose(quat2dcm(β1))*j.jointF[1:3]
    Fb2[4:6] = quat2dcm(β2)*transpose(quat2dcm(β1))*j.jointF[4:6]

    return (Fb1,Fb2)

end
