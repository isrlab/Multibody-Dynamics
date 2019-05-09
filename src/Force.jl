# Function to generate constraint and/or external forces

# include("RigidBody.jl")
# include("OrientationConversion.jl")
# include("Joint.jl")

using LinearAlgebra
using ForwardDiff
using Revise

function ForceCon(j::Joint,extF1::extForces,extF2::extForces,GravityInInertial::Vector{Float64})::Vector{Float64}
    if(j.RB1.m == 0.0)
        # Joint to world. Body is connected to the Inertial frame.
        if j.type == "Revolute"
            Fc = ForceConRevIn(j,extF2,GravityInInertial)
        elseif j.type == "Spherical"
            Fc = ForceConSphIn(j,extF2,GravityInInertial)
        elseif j.type == "Weld"
            Fc = ForceConWeldAllIn(j,extF2,GravityInInertial)
        elseif j.type == "Free"
            Fc = ForceFree(j,extF2,GravityInInertial)
        elseif j.type == "Spring"
            Fc = ForceConSprIn(j,extF2,GravityInInertial)
        else
            error("Joint type not prescribed.")
        end
    else
        if j.type == "Revolute"
            Fc = ForceConRev(j,extF1,extF2,GravityInInertial)
        elseif j.type == "Spherical"
            Fc = ForceConSph(j,extF1,extF2,GravityInInertial)
        elseif j.type == "Weld"
            Fc = ForceConWeldAll(j,extF1,extF2,GravityInInertial)
        elseif j.type == "Spring"
            Fc = ForceConSpr(j,extF1,extF2,GravityInInertial)
        else
            error("Joint type not prescribed.")
        end
    end
    return Fc
end
# For inertial frame and body
function ForceFree(j::Joint,extF::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated due to quaternion constraint
    b2 = j.RB2


    # Revolute Joint has 7 constraints
    A = zeros(1,7); b = [0.0];
    A[1,:] = QuatNormConstraint(j)[1][2,8:14]; b[1] = QuatNormConstraint(j)[2][2]

    # A[1:3,:] = TranslationConstraint(j)[1]; b[1:3] = TranslationConstraint(j)[2]
    # A[5:6,:] = RevJointConstraint(j)[1]; b[5:6] = RevJointConstraint(j)[2]
    # A_in = A[:,8:14]
    # A = zeros(7,14)
    # IMat3 = Matrix{Float64}(I,3,3)
    # t1 = TranslationConstraint(b1.x,rj1)
    # t2 = TranslationConstraint(b2.x,rj2)
    # A[1:3,:] = [-IMat3 t1[1] IMat3 t2[1]]
    # A[4,:] = [zeros(1,11) 1/j.axis[1] -1/j.axis[2] 0]
    # A[5,:] = [zeros(1,11) 1/j.axis[1] 0 -1/j.axis[3]]
    # A[6,:] = [zeros(1,3) transpose(β1) zeros(1,7)]
    # A[7,:] = [zeros(1,10) transpose(β2)]
    #
    # b = zeros(7)
    # b[1:3] = t1[2] + t2[2]
    # b[4:5] = zeros(2)
    # b[6] = -transpose(β1dot)*β1dot
    # b[7] = -transpose(β2dot)*β2dot

    M = genMatM(b2)

    F = genExtF(b2,extF,GravityInInertial)
    Fconstr = ConstraintForceTorque(M,F,A,b)
    return Fconstr
end
function ForceConRevIn(j::Joint, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by revolute joint
    b2 = j.RB2

    # Revolute Joint has 7 constraints
    A = zeros(6,14); b = zeros(6);
    A[1:3,:] = TranslationConstraint(j)[1]; b[1:3] = TranslationConstraint(j)[2]
    A[4,:] = QuatNormConstraint(j)[1][2,:]; b[4] = QuatNormConstraint(j)[2][2]
    A[5:end,:] = RevJointConstraint(j)[1]; b[5:end] = RevJointConstraint(j)[2]

    A_in = A[:,8:14]
    # A = zeros(7,14)
    # IMat3 = Matrix{Float64}(I,3,3)
    # t1 = TranslationConstraint(b1.x,rj1)
    # t2 = TranslationConstraint(b2.x,rj2)
    # A[1:3,:] = [-IMat3 t1[1] IMat3 t2[1]]
    # A[4,:] = [zeros(1,11) 1/j.axis[1] -1/j.axis[2] 0]
    # A[5,:] = [zeros(1,11) 1/j.axis[1] 0 -1/j.axis[3]]
    # A[6,:] = [zeros(1,3) transpose(β1) zeros(1,7)]
    # A[7,:] = [zeros(1,10) transpose(β2)]
    #
    # b = zeros(7)
    # b[1:3] = t1[2] + t2[2]
    # b[4:5] = zeros(2)
    # b[6] = -transpose(β1dot)*β1dot
    # b[7] = -transpose(β2dot)*β2dot

    M = genMatM(b2)

    F = genExtF(b2,extF2,GravityInInertial)

    Fconstr = ConstraintForceTorque(M,F,A_in,b)
    return Fconstr
end

function ForceConWeldAllIn(j::Joint, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by weld joint
    b2 = j.RB2

    # Weld Joint has 8/9 constraints
    A = zeros(9,14); b = zeros(9);
    A[1:3,:] = TranslationConstraint(j)[1]; b[1:3] = TranslationConstraint(j)[2]
    A[4:5,:] = QuatNormConstraint(j)[1]; b[4:5] = QuatNormConstraint(j)[2]
    A[6:9,:] = WeldJointAllConstraint(j)[1]; b[6:9] = WeldJointAllConstraint(j)[2]
    A_in = A[:,8:14]

    M = genMatM(b2)

    F = genExtF(b2,extF2,GravityInInertial)

    Fconstr = ConstraintForceTorque(M,F,A_in,b)
    return Fconstr
end

function ForceConSphIn(j::Joint, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by spherical joint
    b2 = j.RB2

    # Spherical Joint leaves 3 degrees of freedom (remove 4 constraints)
    A = zeros(4,7); b = zeros(4)
    A[1:3,:] = TranslationConstraint(j)[1][:,8:14]; b[1:3] = TranslationConstraint(j)[2]
    A[4,:] = QuatNormConstraint(j)[1][2,8:14]; b[4] = QuatNormConstraint(j)[2][2]

    M = genMatM(b2)

    F = genExtF(b2,extF2,GravityInInertial)

    Fconstr = ConstraintForceTorque(M,F,A,b)
    return Fconstr
end

function ForceConSprIn(j::Joint, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # To compute Spring force and torque
    posSpr1 = j.RB1.x[1:3] + transpose(j.RB1.dcm)*j.pos1 # Inertial position of spring connection on first body
    posSpr2 = j.RB2.x[1:3] + transpose(j.RB2.dcm)*j.pos2 # Inertial position of spring connection on second body

    # First body is the inertial frame

    # Second body
    unitVec = (posSpr2 - posSpr1)/norm(posSpr2-posSpr1,2) # unit Vector in direction of spring (first body to second body)
    F1 = j.k*(norm(posSpr1-posSpr2,2)-j.restLen)*unitVec # Force exerted on first body
    E2 = genE(j.RB2.x[4:7])
    F2 = -F1 # Force exerted by spring on second body
    τ2 = cross(transpose(j.RB2.dcm)*j.pos2,F2) # Torque exerted by spring on second body
    Γb2 = [0.0;j.RB2.dcm*τ2]
    Γu2 = 2*E2*Γb2

    # QuatNormConstraint
    A = zeros(1,7); b = zeros(1)
    A[1,:] = QuatNormConstraint(j)[1][2,8:14]; b[1] = QuatNormConstraint(j)[2][2]

    b2 = j.RB2
    M = genMatM(b2)
    F = genExtF(b2,extF2,GravityInInertial)
    Fconstr = ConstraintForceTorque(M,F,A,b)

    return Fc = [F2+Fconstr[1:3];Γu2+Fconstr[4:7]]
end
# For 2 rigid bodies
function ForceConRev(j::Joint, extF1::extForces, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by revolute joint
    b1 = j.RB1
    b2 = j.RB2

    # Revolute Joint has 7 constraints
    A = zeros(7,14); b = zeros(7);
    A[1:3,:] = TranslationConstraint(j)[1]; b[1:3] = TranslationConstraint(j)[2]
    A[4:5,:] = QuatNormConstraint(j)[1]; b[4:5] = QuatNormConstraint(j)[2]
    A[6:end,:] = RevJointConstraint(j)[1]; b[6:end] = RevJointConstraint(j)[2]
    # A = zeros(7,14)
    # IMat3 = Matrix{Float64}(I,3,3)
    # t1 = TranslationConstraint(b1.x,rj1)
    # t2 = TranslationConstraint(b2.x,rj2)
    # A[1:3,:] = [-IMat3 t1[1] IMat3 t2[1]]
    # A[4,:] = [zeros(1,11) 1/j.axis[1] -1/j.axis[2] 0]
    # A[5,:] = [zeros(1,11) 1/j.axis[1] 0 -1/j.axis[3]]
    # A[6,:] = [zeros(1,3) transpose(β1) zeros(1,7)]
    # A[7,:] = [zeros(1,10) transpose(β2)]
    #
    # b = zeros(7)
    # b[1:3] = t1[2] + t2[2]
    # b[4:5] = zeros(2)
    # b[6] = -transpose(β1dot)*β1dot
    # b[7] = -transpose(β2dot)*β2dot

    M1 = genMatM(b1)
    M2 = genMatM(b2)
    M = [M1 zeros(size(M1))
         zeros(size(M1)) M2]

    F = [genExtF(b1,extF1,GravityInInertial); genExtF(b2,extF2,GravityInInertial)]
    Fconstr = ConstraintForceTorque(M,F,A,b)
    return Fconstr
end

function ForceConWeldAll(j::Joint, extF1::extForces, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by weld joint
    b1 = j.RB1
    b2 = j.RB2

    # Weld Joint has 8/9 constraints
    A = zeros(9,14); b = zeros(9);
    A[1:3,:] = TranslationConstraint(j)[1]; b[1:3] = TranslationConstraint(j)[2]
    A[4:5,:] = QuatNormConstraint(j)[1]; b[4:5] = QuatNormConstraint(j)[2]
    A[6:9,:] = WeldJointAllConstraint(j)[1]; b[6:9] = WeldJointAllConstraint(j)[2]

    M1 = genMatM(b1)
    M2 = genMatM(b2)
    M = [M1 zeros(size(M1))
         zeros(size(M1)) M2]

    F = [genExtF(b1,extF1,GravityInInertial); genExtF(b2,extF2,GravityInInertial)]

    Fconstr = ConstraintForceTorque(M,F,A,b)
    return Fconstr
end

function ForceConSph(j::Joint, extF1::extForces, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by spherical joint
    b1 = j.RB1
    b2 = j.RB2

    # Spherical Joint has 5 constraints
    A = zeros(5,14); b = zeros(5)
    A[1:3,:] = TranslationConstraint(j)[1]; b[1:3] = TranslationConstraint(j)[2]
    A[4:5,:] = QuatNormConstraint(j)[1]; b[4:5] = QuatNormConstraint(j)[2]

    M1 = genMatM(b1)
    M2 = genMatM(b2)
    M = [M1 zeros(size(M1))
         zeros(size(M1)) M2]

    F = [genExtF(b1,extF1,GravityInInertial); genExtF(b2,extF2,GravityInInertial)]

    Fconstr = ConstraintForceTorque(M,F,A,b)
    return Fconstr
end

function ForceConSpr(j::Joint, extF1::extForces, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # To compute Spring force and torque
    posSpr1 = j.RB1.x[1:3] + transpose(j.RB1.dcm)*j.pos1 # Inertial position of spring connection on first body
    posSpr2 = j.RB2.x[1:3] + transpose(j.RB2.dcm)*j.pos2 # Inertial position of spring connection on second body
    unitVec = (posSpr2 - posSpr1)/norm(posSpr2-posSpr1,2) # unit Vector in direction of spring (first body to second body)

    # First body
    E1 = genE(j.RB1.x[4:7])
    F1 = j.k*(norm(posSpr1-posSpr2)-j.restLen)*unitVec # Force exerted by spring on the first body
    τ1 = cross(transpose(j.RB1.dcm)*j.pos1,F1) # Torque exerted by spring on the first body
    Γb1 = [0.0;j.RB1.dcm*τ1]
    Γu1 = 2*E1*Γb1

    # Second body
    E2 = genE(j.RB2.x[4:7])
    F2 = -F1 # Force exerted by spring on second body
    τ2 = cross(transpose(j.RB2.dcm)*j.pos2,F2) # Torque exerted by spring on second body
    Γb2 = [0.0;j.RB2.dcm*τ2]
    Γu2 = 2*E2*Γb2

    # QuatNormConstraint
    A = zeros(2,14); b = zeros(2)
    A = QuatNormConstraint(j)[1]; b = QuatNormConstraint(j)[2]
    b1 = j.RB1
    b2 = j.RB2
    M1 = genMatM(b1)
    M2 = genMatM(b2)
    M = [M1 zeros(size(M1))
         zeros(size(M1)) M2]

    F = [genExtF(b1,extF1,GravityInInertial); genExtF(b2,extF2,GravityInInertial)]
    Fconstr = ConstraintForceTorque(M,F,A,b)

    return Fc = [F1+Fconstr[1:3];Γu1+Fconstr[4:7];F2+Fconstr[8:10];Γu2+Fconstr[11:14]]
end
# forward diff
# Translation constraint
function TranslationConstraint(j::Joint)::Tuple
    # r(y) = transpose(quat2dcm(y))*pos
    # rJac(z) = ForwardDiff.jacobian(r,z)
    # rdot(y::AbstractVector) = ForwardDiff.jacobian(r,y)*βdot

    b1 = j.RB1
    b2 = j.RB2
    rj1 = j.pos1
    rj2 = j.pos2
    β1 = b1.x[4:7]
    β2 = b2.x[4:7]
    β1dot = b1.x[11:14]
    β2dot = b2.x[11:14]

    A = zeros(3,14)
    A[:,1:3] = Matrix{Float64}(I,3,3)
    A[:,8:10] = -Matrix{Float64}(I,3,3)
    A[:,4:7] = TranslationConstraintSupplement(b1.x,rj1)[1]
    A[:,11:14] = -TranslationConstraintSupplement(b2.x,rj2)[1]

    b = zeros(3)

    b[1:3] = -TranslationConstraintSupplement(b1.x,rj1)[2] + TranslationConstraintSupplement(b2.x,rj2)[2]

    return (A,b)
    # rddotRHS = ForwardDiff.jacobian(rdot,β)*βdot
    # rddotLHS = -rJac(β)
    # return (rddotLHS,rddotRHS)
end

function TranslationConstraintSupplement(x::Vector{T},pos::Vector{T}) where T <: Real
    β = x[4:7]
    βdot = x[11:14]
    r(y::Vector{T}) where T<:Real = transpose(quat2dcm(y))*pos
    rJac = z->ForwardDiff.jacobian(r,z)
    rdot(y) = ForwardDiff.jacobian(r,y)*βdot
    rddotRHS = ForwardDiff.jacobian(rdot,β)*βdot
    rddotLHS = rJac(β)
    return (rddotLHS,rddotRHS)
end

function QuatNormConstraint(j::Joint)::Tuple
    A = zeros(2,14)
    b = zeros(2)

    b1 = j.RB1
    b2 = j.RB2

    β1 = b1.x[4:7]
    β2 = b2.x[4:7]
    β1dot = b1.x[11:14]
    β2dot = b2.x[11:14]

    A[1,4:7] = β1
    A[2,11:14] = β2
    b[1] = -transpose(β1dot)*β1dot
    b[2] = -transpose(β2dot)*β2dot
    return (A,b)
end

# function RevJointConstraint(j::Joint)
#     function f(y::Vector{T}) where T <: Real
#         # Computes the Hamilton product of β2 and β1^-1 (i.e. β2*(β1^-1))
#         y1 = [y[1];-y[2:4]]
#         y2 = y[5:8]
#         x = [y2[1]*y1[1] - transpose(y2[2:4])*y1[2:4];
#              y2[1]*y1[2:4] + y1[1]*y2[2:4] + cross(y2[2:4],y1[2:4])]
#         return x
#     end
#     fJac = z->ForwardDiff.jacobian(f,z)
#     fdotJac = z->ForwardDiff.jacobian(fJac,z)
#
#     b1 = j.RB1
#     b2 = j.RB2
#     β1 = b1.x[4:7]
#     β2 = b2.x[4:7]
#     β1dot = b1.x[11:14]
#     β2dot = b2.x[11:14]
#     β = [β1;β2]
#     axis = transpose(quat2dcm(β1))*j.axis
#     a1, a2, a3 = axis
#
#     # β1 = x1[4:7]
#     # β2 = x2[4:7]
#     # β1dot = x1[11:14]
#     # β2dot = x2[11:14]
#
#     y1 = f(β)
#     y2 = fJac(β)
#     y3 = fdotJac(β)
#     x1 = y3[1:16,1:4]*β1dot; x1 = reshape(x1,(4,4))
#     x2 = y3[1:16,5:8]*β2dot; x2 = reshape(x2,(4,4))
#     x3 = y3[17:end,1:4]*β1dot; x3 = reshape(x3,(4,4))
#     x4 = y3[17:end,5:8]*β2dot; x4 = reshape(x4,(4,4))
#     qddotA = [zeros(4,3) -fJac(β)[:,1:4] zeros(4,3) -fJac(β)[:,5:8]]
#     qddotb = (x1+x2)*β1dot + (x3+x4)*β2dot
#
#     nnzAxis = findall(!iszero,axis)
#     zAxis = findall(iszero,axis)
#     l = length(nnzAxis) # Counts the number of nonzero elements in axis
#     A = zeros(2,14); b = zeros(2)
#     if l==3
#         A[1,:] = qddotA[2,:]/a1 - qddotA[3,:]/a2
#         b[1] = qddotb[2]/a1 - qddotb[3]/a2
#         A[2,:] = qddotA[2,:]/a1 - qddotA[4,:]/a3
#         b[2] = qddotb[2]/a1 - qddotb[4]/a3
#     elseif l==2
#         A[1,:] = qddotA[nnzAxis[1]+1,:]/axis[nnzAxis[1]] - qddotA[nnzAxis[2]+1,:]/axis[nnzAxis[2]]
#         b[1] = qddotb[nnzAxis[1]+1]/axis[nnzAxis[1]] - qddotb[nnzAxis[2]+1]/axis[nnzAxis[2]]
#         A[2,:] = qddotA[zAxis[1]+1,:]
#     elseif l==1
#         A[1,:] = qddotA[zAxis[1]+1,:]
#         A[2,:] = qddotA[zAxis[2]+1,:]
#     else
#         error("Joint axis provided not possible.")
#     end
#     return (A,b)
# end

# function RevJointConstraint(j::Joint)
#     # Attempting q_rev(β1,β2) = axis(joint)
#     b1 = j.RB1
#     b2 = j.RB2
#     β1 = b1.x[4:7]
#     β2 = b2.x[4:7]
#     β1dot = b1.x[11:14]
#     β2dot = b2.x[11:14]
#     β = [β1;β2]
#     axis = j.axis
#     function f(y::Vector{T}) where T <: Real
#         # Computes the Hamilton product of β2 and β1^-1 (i.e. β2*(β1^-1))
#         y1 = [y[1];-y[2:4]]
#         y2 = y[5:8]
#         x = [y2[1]*y1[1] - transpose(y2[2:4])*y1[2:4];
#              y2[1]*y1[2:4] + y1[1]*y2[2:4] + cross(y2[2:4],y1[2:4])]
#         return x
#     end
#
#     if(norm(f(β)[2:4],2) > 0)
#         qRev(y) = f(y)[2:4]/norm(f(y)[2:4],2)
#         qRevJac = z->ForwardDiff.jacobian(qRev,z)
#         qRevdotJac = z->ForwardDiff.jacobian(qRevJac,z)
#
#         axisFn(y) = transpose(quat2dcm(y))*axis
#         axisFnJac(y) = ForwardDiff.jacobian(axisFn,y)
#         axisdot(y) = ForwardDiff.jacobian(axisFn,y)*β1dot
#         axisddot = ForwardDiff.jacobian(axisdot,β1)*β1dot
#
#
#         y3 = qRevdotJac(β)
#         x1 = y3[1:12,1:4]*β1dot; x1 = reshape(x1,(3,4))
#         x2 = y3[1:12,5:8]*β2dot; x2 = reshape(x2,(3,4))
#         x3 = y3[13:end,1:4]*β1dot; x3 = reshape(x3,(3,4))
#         x4 = y3[13:end,5:8]*β2dot; x4 = reshape(x4,(3,4))
#
#
#         A = zeros(3,14); b = zeros(3)
#         A[:,4:7] = -axisFnJac(β1) + qRevJac(β)[:,1:4]
#         A[:,11:14] = qRevJac(β)[:,5:8]
#         b = axisddot - (x1+x2)*β1dot - (x3+x4)*β2dot
#     else
#         A = zeros(3,14); b = zeros(3)
#     end
#
#     # y1 = f(β)
#     # y2 = fJac(β)
#     # y3 = fdotJac(β)
#     # x1 = y3[1:16,1:4]*β1dot; x1 = reshape(x1,(4,4))
#     # x2 = y3[1:16,5:8]*β2dot; x2 = reshape(x2,(4,4))
#     # x3 = y3[17:end,1:4]*β1dot; x3 = reshape(x3,(4,4))
#     # x4 = y3[17:end,5:8]*β2dot; x4 = reshape(x4,(4,4))
#     # qddotA = [zeros(4,3) -fJac(β)[:,1:4] zeros(4,3) -fJac(β)[:,5:8]]
#     # qddotb = (x1+x2)*β1dot + (x3+x4)*β2dot
#     #
#     # nnzAxis = findall(!iszero,axis)
#     # zAxis = findall(iszero,axis)
#     # l = length(nnzAxis) # Counts the number of nonzero elements in axis
#     # A = zeros(2,14); b = zeros(2)
#     # if l==3
#     #     A[1,:] = qddotA[2,:]/a1 - qddotA[3,:]/a2
#     #     b[1] = qddotb[2]/a1 - qddotb[3]/a2
#     #     A[2,:] = qddotA[2,:]/a1 - qddotA[4,:]/a3
#     #     b[2] = qddotb[2]/a1 - qddotb[4]/a3
#     # elseif l==2
#     #     A[1,:] = qddotA[nnzAxis[1]+1,:]/axis[nnzAxis[1]] - qddotA[nnzAxis[2]+1,:]/axis[nnzAxis[2]]
#     #     b[1] = qddotb[nnzAxis[1]+1]/axis[nnzAxis[1]] - qddotb[nnzAxis[2]+1]/axis[nnzAxis[2]]
#     #     A[2,:] = qddotA[zAxis[1]+1,:]
#     # elseif l==1
#     #     A[1,:] = qddotA[zAxis[1]+1,:]
#     #     A[2,:] = qddotA[zAxis[2]+1,:]
#     # else
#     #     error("Joint axis provided not possible.")
#     # end
#     return (A,b)
# end

# function RevJointConstraint(j::Joint)
#     # Attempting βr = β2*β1Inv; equating axes in inertial frame
#     b1 = j.RB1
#     b2 = j.RB2
#     β1 = b1.x[4:7]
#     β2 = b2.x[4:7]
#     β1dot = b1.x[11:14]
#     β2dot = b2.x[11:14]
#     β = [β1;β2]
#     axis = j.axis
#     function f(y::Vector{T}) where T <: Real
#         # Returns the Euler axis in the body frame of b1 of the revolute joint rotation quaternion
#         y1 = y[1:4]
#         y1inv = [y[1];-y[2:4]]
#         y2 = y[5:8]
#         # β2inβ1inv= [y2[1];quat2dcm(β1inv)*transpose(quat2dcm(y2))*y2[2:4]]
#         yRev = quaternionProduct(y2,y1inv)
#         axisRev = transpose(quat2dcm(y1))*yRev[2:4]/norm(yRev[2:4],2)
#         # βRev = dcm2quat(CRev)
#         # axisRev = transpose(quat2dcm(y1))*βRev[2:4]/norm(βRev[2:4],2)
#         return axisRev
#     end
#     fJac = z->ForwardDiff.jacobian(f,z)
#     fdotJac = z->ForwardDiff.jacobian(fJac,z)
#
#     axisFn(y) = transpose(quat2dcm(y))*axis
#     axisFnJac(y) = ForwardDiff.jacobian(axisFn,y)
#     axisdot(y) = ForwardDiff.jacobian(axisFn,y)*β1dot
#
#     β1inv = [β[1];-β[2:4]]
#     βRev = quaternionProduct(β2,β1inv)
#     if norm(βRev[2:4],2) > 0
#         y3 = fdotJac(β)
#         x1 = y3[1:12,1:4]*β1dot; x1 = reshape(x1,(3,4))
#         x2 = y3[1:12,5:8]*β2dot; x2 = reshape(x2,(3,4))
#         x3 = y3[13:24,1:4]*β1dot; x3 = reshape(x3,(3,4))
#         x4 = y3[13:24,5:8]*β2dot; x4 = reshape(x4,(3,4))
#
#         A = zeros(3,14); b = zeros(3)
#         axisddot = ForwardDiff.jacobian(axisdot,β1)*β1dot
#         A[:,4:7] = fJac(β)[:,1:4]  - axisFnJac(β1)
#         A[:,11:14] = fJac(β)[:,5:8]
#         b = - (x1+x2)*β1dot - (x3+x4)*β2dot + axisddot
#     else
#         A = zeros(3,14); b = zeros(3)
#     end
#     return (A,b)
# end

# function RevJointConstraint(j::Joint)
#     # Attempting from chapter on Spatial Kinematics
#     b1 = j.RB1
#     b2 = j.RB2
#     β1 = b1.x[4:7]
#     β2 = b2.x[4:7]
#     β1dot = b1.x[11:14]
#     β2dot = b2.x[11:14]
#     axis = j.axis
#     rj2 = j.pos2
#     βaug = [β1;β2;axis] # β augmented.
#     function f(y::Vector{T}) where T <: Real
#         # y = [β1;β2;axis]
#         # Returns the cross product of the joint axis expressed in the frames of the two bodies
#         y1 = y[1:4]
#         y2 = y[5:8]
#         # rj2 = y[9:11]
#         v1 = transpose(quat2dcm(y1))*axis
#         v2 = transpose(quat2dcm(y2))*axis
#         return skewX(v1)*v2
#     end
#     fJac = z->ForwardDiff.jacobian(f,z)
#     fdotJac = z->ForwardDiff.jacobian(fJac,z)
#
#     y3 = fdotJac(βaug)
#     x1 = y3[1:12,1:4]*β1dot; x1 = reshape(x1,(3,4))
#     x2 = y3[1:12,5:8]*β2dot; x2 = reshape(x2,(3,4))
#     x3 = y3[13:24,1:4]*β1dot; x3 = reshape(x3,(3,4))
#     x4 = y3[13:24,5:8]*β2dot; x4 = reshape(x4,(3,4))
#
#     A = zeros(3,14); b = zeros(3)
#     A[:,4:7] = fJac(βaug)[:,1:4]
#     A[:,11:14] = fJac(βaug)[:,5:8]
#     b = - (x1+x2)*β1dot - (x3+x4)*β2dot
#
#     return (A,b)
# end

# function RevJointConstraint(j::Joint)
#     # Attempting β2axis = jt.axis in Inertial
#     b1 = j.RB1
#     b2 = j.RB2
#     β1 = b1.x[4:7]
#     β2 = b2.x[4:7]
#     β1dot = b1.x[11:14]
#     β2dot = b2.x[11:14]
#     axis = j.axis
#     β1aug = [β1;axis] # β augmented.
#
#     if norm(β2[2:4],2) > 0
#         function f(β2::Vector{T}) where T <: Real
#             # Returns the Euler axis of β2
#             return β2[2:4]/norm(β2[2:4],2)
#         end
#         fJac = z->ForwardDiff.jacobian(f,z)
#         f2 = z->fJac(z)*β2dot
#         f2Jac = z->ForwardDiff.jacobian(f2,z)
#
#         function g(β1Aug::Vector{T}) where T <: Real
#             # Returns the joint axis in the inertial frame
#             return transpose(quat2dcm(β1Aug[1:4]))*β1Aug[5:7]
#         end
#         gJac= z->ForwardDiff.jacobian(g,z)
#         g2 = z->gJac(z)[:,1:4]*β1dot
#         g2Jac = z->ForwardDiff.jacobian(g2,z)
#
#         A = zeros(3,14); b = zeros(3)
#         A[:,4:7] = -gJac(β1aug)[:,1:4]
#         A[:,11:14] = fJac(β2)
#         b = -f2Jac(β2)*β2dot + g2Jac(β1aug)[:,1:4]*β1dot
#     else
#         A = zeros(3,14); b = zeros(3)
#     end
#     @show A
#     @show b
#     return (A,b)
# end

# function RevJointConstraint(j::Joint)
#     # Attempting βRevaxis = axis(dcm2quat(C(β2)*R(β1)))
#     b1 = j.RB1
#     b2 = j.RB2
#     β1 = b1.x[4:7]
#     β2 = b2.x[4:7]
#     β1dot = b1.x[11:14]
#     β2dot = b2.x[11:14]
#     axis = j.axis
#     β = [β1;β2] # β augmented.
#
#     function f(y::Vector{T}) where T <: Real
#         y1 = y[1:4]; y2 = y[5:8]
#         CβRev = quat2dcm(y2)*transpose(quat2dcm(y1))
#         βRev = dcm2quat(CβRev)
#         ax = βRev[2:4]/norm(βRev[2:4],2)
#         return ax
#     end
#     fJac = z->ForwardDiff.jacobian(f,z)
#     fdotJac = z->ForwardDiff.jacobian(fJac,z)
#
#     y3 = fdotJac(β)
#     x1 = y3[1:12,1:4]*β1dot; x1 = reshape(x1,(3,4))
#     x2 = y3[1:12,5:8]*β2dot; x2 = reshape(x2,(3,4))
#     x3 = y3[13:24,1:4]*β1dot; x3 = reshape(x3,(3,4))
#     x4 = y3[13:24,5:8]*β2dot; x4 = reshape(x4,(3,4))
#
#     A = zeros(3,14); b = zeros(3)
#     @show f(β)
#     if !any(i->isnan(f(β)[i]),1:3) # making sure there are no NaN values in the vector
#         A[:,4:7] = fJac(β)[:,1:4]
#         A[:,11:14] = fJac(β)[:,5:8]
#         b = - (x1+x2)*β1dot - (x3+x4)*β2dot
#     end
#     @show A
#     @show b
#     return (A,b)
# end

function RevJointConstraint(j::Joint)::Tuple
    # Attempting ω of body 2 in body 1 frame only has component about the joint axis
    b1 = j.RB1
    b2 = j.RB2
    β1 = b1.x[4:7]
    β2 = b2.x[4:7]
    β1dot = b1.x[11:14]
    β2dot = b2.x[11:14]
    axis = j.axis
    βaug = [β1;β2;β2dot]
    function f(y::Vector{T}) where T <: Real
        y1  = y[1:4]; #β1
        y2 = y[5:8]; y2dot = y[9:12] #β2;β2dot
        ωb2 = angVel(y2,y2dot)
        ωb2inb1 = quat2dcm(y1)*transpose(quat2dcm(y2))*ωb2
        return ωb2inb1
    end
    fJac = z->ForwardDiff.jacobian(f,z)

    Ax = zeros(3,14); bx = zeros(3)
    Ax[:,11:14] = fJac(βaug)[:,9:12]
    bx = -fJac(βaug)[:,1:4]*β1dot - fJac(βaug)[:,5:8]*β2dot

    zAxis = findall(iszero,axis)
    A = zeros(2,14); b = zeros(2)
    A[1,:] = Ax[zAxis[1],:]
    A[2,:] = Ax[zAxis[2],:]
    b[1] = bx[zAxis[1]]
    b[2] = bx[zAxis[2]]

    return (A,b)

end

function WeldJointAllConstraint(j::Joint)::Tuple
    b1 = j.RB1
    b2 = j.RB2

    A = zeros(4,14)
    b = zeros(4)
    y = WeldJointAllConstraintSupplement(b1.x,b2.x)
    A[:,4:7] = y[1][1]
    A[:,11:14] = y[1][2]
    b = y[2]
    return (A,b)
end

function WeldJointAllConstraintSupplement(x1::Vector{T},x2::Vector{T}) where T <: Real
    function f(y::Vector{T}) where T <: Real
        # Computes the Hamilton product of β2 and β1^-1 (i.e. β2*()β1^-1))
        y1 = [y[1];-y[2:4]]
        y2 = y[5:8]
        x = [y2[1]*y1[1] - transpose(y2[2:4])*y1[2:4];
             y2[1]*y1[2:4] + y1[1]*y2[2:4] + cross(y2[2:4],y1[2:4])]
        return x
    end
    fJac = z->ForwardDiff.jacobian(f,z)
    fdotJac = z->ForwardDiff.jacobian(fJac,z)

    β1 = x1[4:7]
    β2 = x2[4:7]
    β1dot = x1[11:14]
    β2dot = x2[11:14]
    β = [β1;β2]

    y1 = f(β)
    y2 = fJac(β)
    y3 = fdotJac(β)
    y3r = reshape(y3,(4,8,8))
    x1 = y3[1:16,1:4]*β1dot; x1 = reshape(x1,(4,4))
    x2 = y3[1:16,5:8]*β2dot; x2 = reshape(x2,(4,4))
    x3 = y3[17:end,1:4]*β1dot; x3 = reshape(x3,(4,4))
    x4 = y3[17:end,5:8]*β2dot; x4 = reshape(x4,(4,4))
    # β = rand(8); β1dot = rand(4); β2dot = rand(4);
    LHS = (-fJac(β)[:,1:4],-fJac(β)[:,5:8])
    RHS = (x1+x2)*β1dot + (x3+x4)*β2dot
    return (LHS,RHS)

end
function ConstraintForceTorque(M::Matrix{Float64},F::Vector{Float64},A::Matrix{Float64},b::Vector{Float64})::Vector{Float64}
    # @show A
    # @show b
    # @show A*inv(M)*F
    if !isreal(sqrt(M)) # insignificant imaginary values showing up
        @show isreal(sqrt(M))
        @show sqrt(M)
        M1 = real(sqrt(M))
        M2 = real(M^(-0.5))
        Fc = M1*pinv(A*M2)*(b-A*inv(M)*F)
    else
        Fc = sqrt(M)*pinv(A*M^(-0.5))*(b-A*inv(M)*F)
    end
    return Fc
    # Fc is in the inertial frame
end

function genMatM(b::RigidBody)::Matrix{Float64}
    # Function to generate mass matrix for each rigid body
    β = b.x[4:7]
    E = genE(β)
    J = b.J
    M = [b.m*Matrix{Float64}(I,3,3)          zeros(3,4)
                         zeros(4,3) 4*transpose(E)*J*E]
end

function genExtF(b::RigidBody,extF::extForces,GravityInInertial::Vector{Float64})::Vector{Float64}
    # Function to generate augmented external Force vector for unconstrained system
    # External Forces are always in the body frame
    β = b.x[4:7]
    βdot = b.x[11:14]
    b.dcm = quat2dcm(β)
    b.ω = angVel(β,βdot)
    E = genE(β)
    Edot = -genE(βdot)
    TotalMoment = zeros(3)
    # J0 = rand(1)# random positive number
    for i in 1:size(extF.Forces)[1]
        TotalMoment = TotalMoment + cross(extF.Positions[i,:],extF.Forces[i,:])
    end
    TotalMoment = TotalMoment + sum(extF.Torques,dims=1)[:]
    # TotalMoment = transpose(b.dcm)*TotalMoment
    Γb = [0.0;TotalMoment] # In the body frame
    Γu = 2*E*Γb
    # @show β
    # @show transpose(b.dcm)
    # @show (sum(extF.Forces,dims=1)[:])

    F = [transpose(b.dcm)*(sum(extF.Forces,dims=1)[:]) + b.m*GravityInInertial
         Γu - 8*transpose(Edot)*b.J*E*βdot - 4*b.J[1,1]*(transpose(βdot)*βdot)*β]
    return F
    # F returned is in the inertial frame
end