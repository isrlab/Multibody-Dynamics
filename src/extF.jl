#
include("RigidBody.jl")
include("Joint.jl")
include("OrientationConversion.jl")

function extF(t::T, j::Vector{Joint})::Vector{extForces} where T<:Real
    # Function to generate external forces
    # that may or may not be functions of time
    extFList = zeroExtForceVec(length(j)+1)

    # For reference, calculate total mass of robot
    TotalMass = (0.0);
    for i=1:length(j)
        TotalMass += j[i].RB2.m
    end

    # Pendulum
    # extFList[2] = extForces([0.0 0 0.1],[0.0 0 0],[0.0 0.0 0.0]);

    # Double Pendulum Test
    # extFList[2] = zeroExtForce()
    # extFList[3] = extForces([0.0 0 0],[0.0 0 0],[0.0 0.0 0.0])

    # Cube Prop Test
    # extFList[2] = extForces([0.0 0.0 0.0], zeros(1,3), [0.0 0.0 0.0])
    # extFList[3] = extForces([0.0 0.0 0.0], zeros(1,3), [0.0 0.0 0.1])

    # QuadTest
    # extFList[3] = extForces([0.0 0.0 0.0], zeros(1,3), [0.0 0.0 0.1])#*sin(t)])
    # extFList[4] = extForces([0.0 0.0 0.0], zeros(1,3), [0.0 0.0 0.1])#*sin(t)])
    # extFList[5] = extForces([0.0 0.0 0.0], zeros(1,3), [0.0 0.0 -0.1*sin(t)])
    # extFList[6] = extForces([0.0 0.0 0.0], zeros(1,3), [0.0 0.0 -0.1*sin(t)])

    nRotors = length(j) - 1;
    C_T = 1.27e-6;
    Ω_h = sqrt(TotalMass*(9.81)/nRotors/C_T);
    flt = ones(nRotors);  # fault vector for the n Rotors
    # flt = [1.0;0.5;1.0;0.5];
    for i=3:j[end].RB2.bodyID
        extFList[i] = genRotorF(Ω_h,flt[i-2]);
    end

    # grav = [0.0,0.0,-9.806]

    return extFList
end

function genJointF(t::T,j::Joint)::Tuple{Vector{T}, Vector{T}} where T<:Real
    # For actuated joints.
    # j.jointF gives us the actuation force and torque applied ...
    # ... on the child link in the body frame of the parent link
    β1 = j.RB1.x[4:7]
    β2 = j.RB2.x[4:7]

    Fb1 = -j.jointF # Joint Force acting on parent link
    Fb2 = Vector{T}(undef,6) # Force acting on child link
    Fb2[1:3] = quat2dcm(β2)*transpose(quat2dcm(β1))*j.jointF[1:3]
    Fb2[4:6] = quat2dcm(β2)*transpose(quat2dcm(β1))*j.jointF[4:6]

    if j==2
        print("Fb = \t")
        println(Fb2)
    end
    return (Fb1,Fb2)

end

function genRotorF(Ω,f)
    # f: fault signal ∈ [0,1]. 0 indicates total fault and 1 indicates nominal performance
    C_T = 1.27e-6; # thrust coefficient
    T = C_T*Ω^2; # thrust
    F = f*[0;0;T];
    prop_Fτ = extForces(permutedims(F), zeros(1,3), zeros(1,3));
    return prop_Fτ;
end

## old stuff
# function genMotorF(t::Float64,j::Joint)::Tuple{Vector{Float64}, Vector{Float64}}
#     # For prop joints.
#     # input PWM value from the user (to be supplied from a different function in the future)
#     # lookupTable fn gives us the actuation force and torque applied
#     # ... on the child link in its body frame
#     β1 = j.RB1.x[4:7]
#     β2 = j.RB2.x[4:7]
#
#     FArr =     [1000      0      0;
#                 1100      0      0;
#                 1200 0.0261 0.1564;
#                 1300 0.0719 0.4443;
#                 1400 0.1147 0.7171;
#                 1500 0.1585 0.9877;
#                 1600 0.2138 1.3191]
#
#     pwm = 1400
#     # println("Enter PWM.")
#     # inp = parse(Int64,readline())
#     # if (inp==0) pwm = 1000 else pwm = inp end
#
#     Fb1 = zeros(6) # Joint Force acting on parent link
#     Fb2 = zeros(6) # Joint Force acting on child link
#
#     # Get Force and torque generated on the props in the prop frame
#     Fb2[3], Fb2[6] = lookupTable(FArr,pwm)
#     Fb1[1:3] = quat2dcm(β1)*transpose(quat2dcm(β2))*(-Fb2)[1:3]
#     Fb1[4:6] = quat2dcm(β1)*transpose(quat2dcm(β2))*(-Fb2)[4:6]
#
#     return (Fb1,Fb2)
# end
# ##
# # Testing lookup table implementation for PWM to F,τ
#
# # const τArr = [1000      0;
# #         1100      0;
# #         1200 0.1564;
# #         1300 0.4443;
# #         1400 0.7171;
# #         1500 0.9877;
# #         1600 1.3191]
#
# function lookupTable(FArr::Matrix{Float64},pwm::Real):: Tuple{Float64,Float64}
#
#     a = findfirst(isequal(pwm),FArr)
#     F = FArr[a[1],2]
#     τ = FArr[a[1],3]
#
#     return (F,τ)
# end
