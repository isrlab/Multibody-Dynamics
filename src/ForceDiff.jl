# Function to generate constraint and/or external forces

# include("RigidBody.jl")
# include("OrientationConversion.jl")
# include("Joint.jl")
using GenericLinearAlgebra
using LinearAlgebra
using ForwardDiff
using Revise
using StaticArrays

# function Constraint(x::AbstractArray{T}, j::Tuple{Vararg{Joint}}, extFList::Vector{extForces}, GravityInInertial::MArray{Tuple{3},Real,1,3})::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
#
#     Fc, Ffinal = FcFn(x, j, extFList, GravityInInertial);
#
#     # jb = ForwardDiff.jacobian(z -> FcFn([x[1:14];z], j, extFList, GravityInInertial)[1],x[15:end])
#     # println("jb(Fc) = ", jb)
#     # sleep(1000);
#
#     # Partition generated ForceConstr for each rigid body.
#     ForceConstr = zeros(T,(7,length(j)+1))
#     for i=1:length(j)
#         k = j[i].RB2.bodyID
#         ForceConstr[:,k] = Fc[7*(k-2)+1:7*(k-1)]
#     end
#
#     # Partition generated Ffinal (unconstrainedF) for each rigid body.
#     unconstrF = zeros(T,(7,length(j)+1))
#     for i=1:length(j)
#         k = j[i].RB2.bodyID
#         unconstrF[:,k] = Ffinal[7*(k-2)+1:7*(k-1)]
#     end
#     # Return unconstrF, ForceConstr.
#     return (unconstrF,ForceConstr)
#
# end
#
# function FcFn(x::AbstractArray{T}, j::Tuple{Vararg{Joint}}, extFList::Vector{extForces}, GravityInInertial::MArray{Tuple{3},Real,1,3}) where T<:Real
#
#     Afinal, bfinal = Ab_VecOfMat(x,j);
#
#     ## Check Afinal jac
#     # jb = ForwardDiff.jacobian(z -> Ab_VecOfMat(z,j)[1],x);
#     # println("jb = ", jb[:,15:end])
#     # sleep(1000);
#
#     ## Check bfinal jac
#     # jb = ForwardDiff.jacobian(z -> Ab_VecOfMat(z,j)[2],x);
#     # println("jb(bfinal) = ", jb[:,15:28])
#     # sleep(1000);
#
#     Ffinal = assembleF(x,j,extFList,GravityInInertial)
#     ## Check Ffinal jac
#     # jb = ForwardDiff.jacobian(z -> assembleF(z,j,extFList, GravityInInertial),x)
#     # println("jb = ", jb[:,15:28])
#     # sleep(1000);
#
#     Mfinal, Ms, Ms_inv = M_mat(x,j);
#     ## Check Mfinal jac
#     # jb = ForwardDiff.jacobian(z -> M_mat(z,j)[3], x)
#     # println("jb(Ms_inv) = ", jb[:,15:28])
#     # sleep(1000);
#
#     accUnconstr = Mfinal\Ffinal # Unconstrained Acceleration
#
#     if !isreal(Ms) # insignificant imaginary values showing up
#         M1 = real(Ms)
#         # println("Imag part magnitude =", norm(imag(Ms)))
#         M2 = real(Ms_inv)
#
#         tempMat = Afinal*M2
#         tempMat_pinv = tempMat'*((tempMat*tempMat')\I(size(tempMat)[1]));
#         Fc = Ms*tempMat_pinv*(bfinal - Afinal*accUnconstr);
#     else
#         tempMat = Afinal*Ms_inv;
#         tempMat_pinv = tempMat'*((tempMat*tempMat')\I(size(tempMat)[1]));
#         Fc = Ms*tempMat_pinv*(bfinal - Afinal*accUnconstr); # using Udwadia's formulation
#     end
#     # Fc and Ffinal are in the inertial frame
# return Fc, Ffinal
# end

function Constraint(x::AbstractArray{T}, u::AbstractArray{S}, j::Vector{Joint},  GravityInInertial::Vector{Float64})::Tuple{AbstractArray{Union{T,S}}, AbstractArray{Union{T,S}}} where {T<:Real, S<:Real}

    Fc, Ffinal = FcFn(x, u, j, GravityInInertial);

    # jb = ForwardDiff.jacobian(z -> FcFn([x[1:14];z], j, extFList, GravityInInertial)[1],x[15:end])
    # println("jb(Fc) = ", jb)
    # sleep(1000);

    # Partition generated ForceConstr for each rigid body.
    ForceConstr = Matrix{Union{T,S}}(undef,(7,length(j)+1))
    ForceConstr[:,1] = zeros(7);
    # ForceConstr = zeros(S,(7,length(j)+1))
    for i=1:length(j)
        k = j[i].RB2.bodyID
        ForceConstr[:,k] = Fc[7*(k-2)+1:7*(k-1)]
    end

    # Partition generated Ffinal (unconstrainedF) for each rigid body.
    unconstrF = Matrix{Union{T,S}}(undef, (7,length(j)+1));
    unconstrF[:,1] = zeros(7);
    # unconstrF = zeros(S,(7,length(j)+1))
    for i=1:length(j)
        k = j[i].RB2.bodyID
        unconstrF[:,k] = Ffinal[7*(k-2)+1:7*(k-1)]
    end
    # Return unconstrF, ForceConstr.
    return (unconstrF,ForceConstr)

end

function FcFn(x::AbstractArray{T}, u::AbstractArray{S}, j::Vector{Joint}, GravityInInertial::Vector{Float64}) where {T<:Real, S<:Real}

    Afinal, bfinal = Ab_VecOfMat(x,j);

    ## Check Afinal jac
    # jb = ForwardDiff.jacobian(z -> Ab_VecOfMat(z,j)[1],x);
    # println("jb = ", jb[:,15:end])
    # sleep(1000);

    ## Check bfinal jac
    # jb = ForwardDiff.jacobian(z -> Ab_VecOfMat(z,j)[2],x);
    # println("jb(bfinal) = ", jb[:,15:28])
    # sleep(1000);

    Ffinal = assembleF(x,u,j,GravityInInertial)
    ## Check Ffinal jac
    # jb = ForwardDiff.jacobian(z -> assembleF(z,j,extFList, GravityInInertial),x)
    # println("jb = ", jb[:,15:28])
    # sleep(1000);

    Mfinal, Ms, Ms_inv = M_mat(x,j);
    ## Check Mfinal jac
    # jb = ForwardDiff.jacobian(z -> M_mat(z,j)[3], x)
    # println("jb(Ms_inv) = ", jb[:,15:28])
    # sleep(1000);

    accUnconstr = Mfinal\Ffinal # Unconstrained Acceleration

    if !isreal(Ms) # insignificant imaginary values showing up
        M1 = real(Ms)
        # println("Imag part magnitude =", norm(imag(Ms)))
        M2 = real(Ms_inv)

        tempMat = Afinal*M2
        tempMat_pinv = tempMat'*((tempMat*tempMat')\I(size(tempMat)[1]));
        Fc = Ms*tempMat_pinv*(bfinal - Afinal*accUnconstr);
    else
        tempMat = Afinal*Ms_inv;
        tempMat_pinv = tempMat'*((tempMat*tempMat')\I(size(tempMat)[1]));
        Fc = Ms*tempMat_pinv*(bfinal - Afinal*accUnconstr); # using Udwadia's formulation
    end
    # Fc and Ffinal are in the inertial frame
return Fc, Ffinal
end


function Ab_VecOfMat(x::AbstractArray{T}, j::Vector{Joint})::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    nCols = 7*(j[end].RB2.bodyID-1);
    nConstrEqns_j = zeros(Integer,length(j));
    ## Number of constraint equations for each joint in robot
    for i=1:length(j)
        # Use jointDict to identify number of nConstrEqns in each joint
        if j[i].RB1.bodyID > 1 # First body not the inertial frame
            nConstrEqns_j[i] = get(JointDict, j[i].type, 0)
        else
            nConstrEqns_j[i]= get(JointDictIn, j[i].type, 0)
        end
    end
    nConstrEqns = sum(nConstrEqns_j) # Total number of constraint equations in

    Afinal = zeros(T,(nConstrEqns,nCols))
    bfinal = zeros(T,(nConstrEqns))

    for i = 1:length(j)
        if i > 1
            b1 = j[i].RB1.bodyID-1; b2 = j[i].RB2.bodyID-1;
            nR = sum(nConstrEqns_j[1:i-1])
            Afinal[nR+1:nR + nConstrEqns_j[i],7*(b1-1)+1:7*b1] = ForceCon(x,j[i])[1][:,1:7] # not required for first joint
            Afinal[nR+1:nR + nConstrEqns_j[i],7*(b2-1)+1:7*b2] = ForceCon(x,j[i])[1][:,8:14]
            bfinal[nR+1:nR + nConstrEqns_j[i]] = ForceCon(x,j[i])[2]

            ## old
            # Afinal[nConstrEqns_j[i-1]+1:nConstrEqns_j[i-1] + nConstrEqns_j[i],7*(b1-1)+1:7*b1] = ForceCon(x,j[i])[1][:,1:7] # not required for first joint
            # Afinal[nConstrEqns_j[i-1]+1:nConstrEqns_j[i-1] + nConstrEqns_j[i],7*(b2-1)+1:7*b2] = ForceCon(x,j[i])[1][:,8:14]
            # bfinal[nConstrEqns_j[i-1]+1:nConstrEqns_j[i-1] + nConstrEqns_j[i]] = ForceCon(x,j[i])[2]
        else
            b2 = j[i].RB2.bodyID-1;
            Afinal[1:nConstrEqns_j[i],7*(b2-1)+1:7*b2] = ForceCon(x,j[i])[1][:,8:14]
            bfinal[1:nConstrEqns_j[i]] = ForceCon(x,j[i])[2]
        end

    end

    # Quaternion constraint gets repeated for bodies common to multiple joints
    Afinal, ids = unique_ids(Afinal);
    bfinal = bfinal[ids];
    return Afinal, bfinal
end

function M_mat(x::AbstractArray{T}, j::Vector{Joint})::Tuple{AbstractArray{T}, AbstractArray{T}, AbstractArray{T}} where T<:Real
    Mfinal = assembleM(x,j);
    if typeof(x[1]) == Float64 # eigendecomposition not required for sqrt
        Ms = real(sqrt(Mfinal)) # insignificant imaginary values showing up
        Ms_inv = Ms\I(size(Ms)[1])
    else # Dual type variable (for ForwardDiff), applicable when taking jacobian
        # println("type of Mfinal = ", typeof(Mfinal))
        vals, vecs = eigen(Mfinal);
        # println("typeof(vals) = ", typeof(vals))
        vals_sq = zeros(T,size(vals))
        # vals_sq_inv = zeros(T,size(vals))
        for i=1:length(vals)
            vals_sq[i] = sqrt(vals[i])
        end
        Ms = vecs*diagm(vals_sq)*(vecs\I(size(vecs,1))) #sqrt(M)
        Ms_inv = Ms\I(size(Ms)[1]); #M^(-0.5)
    end
    return Mfinal, Ms, Ms_inv
end

function assembleM(x::AbstractArray{T}, j::Vector{Joint})::AbstractArray{T} where T<:Real
    maxBodyID = j[end].RB2.bodyID
    v = collect(3:maxBodyID)

    Mfinal = zeros(T,(7*(maxBodyID-1),7*(maxBodyID-1)))

    # First body always connected to inertial frame
    Mfinal[1:7,1:7] = genMatM(x,j[1].RB2)

    # Iterate over joints
    for i = 2:length(j)
        b1 = j[i].RB1; b2 = j[i].RB2;
        if b2.bodyID in v
            Mfinal[7*(i-1)+1:7*i,7*(i-1)+1:7*i] = genMatM(x,b2)
            # Remove corresponding element from v
            filter!(x->x≠b2.bodyID,v)
        end
    end
    return Mfinal
end

# function assembleF(x::AbstractArray{T}, j::Tuple{Vararg{Joint}}, extFList::Vector{extForces}, GravityInInertial::MArray{Tuple{3},Real,1,3})::AbstractArray{T} where T<:Real
#     maxBodyID = j[end].RB2.bodyID
#     F = AbstractArray{T}(undef,7*(maxBodyID-1))
#
#     for i=1:length(j)
#         k = j[i].RB2.bodyID - 1
#         F[7*(k-1)+1:7*k] =
#         genExtF(x,j[i].RB2,extFList[k+1],GravityInInertial)
#     end
#     # jb = ForwardDiff.jacobian(z -> genExtF(z,j[1].RB2,extFList[2], GravityInInertial),x)
#     # println("jb = ", jb[:,15:28])
#     # sleep(1000);
#     return F
# end

function assembleF(x::AbstractArray{T}, u::AbstractArray{S}, j::Vector{Joint}, GravityInInertial:: Vector{Float64})::AbstractArray{Union{T,S}} where {T<:Real, S<:Real}
    maxBodyID = j[end].RB2.bodyID
    F = Vector{Union{T,S}}(undef,7*(maxBodyID-1))

    for i=1:length(j)
        k = j[i].RB2.bodyID - 1
        F[7*(k-1)+1:7*k] =
        genExtF(x,j[i].RB2,u[:,k+1],GravityInInertial)

        if j[i].type == "Spring"
            F_b1, F_b2 = ForceConSpr(x,j[i]); # Forces acting on the bodies connected by the spring. Note that a spring does not constrain 2 bodies any more than they previously were, i.e., degrees of freedom remain the same.

            b1_id = j[i].RB1.bodyID -1; b2_id = j[i].RB2.bodyID - 1;
            if j[i].RB1.m == 0.0 # First body is the inertial frame
                F[7*(b2_id-1)+1:7*b2_id] += F_b2;
            else # First body is not the inertial frame
                F[7*(b1_id-1)+1:7*b1_id] += F_b1;
                F[7*(b2_id-1)+1:7*b2_id] += F_b2;
            end
        end
    end
    return F
end

function ForceCon(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    if(j.RB1.m == 0.0)
        # Joint to world. Body is connected to the Inertial frame.
        if j.type == "Revolute"
            A,b = ForceConRevIn(x,j)
            # function tempFn(x::AbstractArray{T}, j::Joint) where T<:Real
            #     Ar = ForceConRevIn(x,j)[1]
            #     Arp = pinv(Ar)
            #     return Arp
            # end
            # jb = ForwardDiff.jacobian(z -> tempFn([x[1:14];z],j),x[15:28]);
            # println("jb(Ar) = ", jb)#[:,15:28])
            # sleep(1000);
        elseif j.type == "Revolute2"
            A,b = ForceConRev2In(x,j)
        elseif j.type == "Spherical"
            A,b = ForceConSphIn(x,j)
            # function tempFn2(x::AbstractArray{T}, j::Joint) where T<:Real
            #     Ar = ForceConSphIn(x,j)[1]
            #     Arp = pinv(Ar)
            #     return Arp
            # end
            # jb = ForwardDiff.jacobian(z -> tempFn2([x[1:14];z],j),x[15:28]);
            # println("jb(Ar) = ", jb)#[:,15:28])
            # sleep(1000);
        elseif j.type == "Weld"
            A,b = ForceConWeldAllIn(x,j)
        elseif j.type == "Free" || j.type == "Spring"
            A,b = ForceFreeIn(x,j)
        else
            error("Joint type not prescribed.")
        end
    else
        if j.type == "Revolute"
            A,b = ForceConRev(x,j)
        elseif j.type == "Revolute2"
            A,b = ForceConRev2(x,j)
        elseif j.type == "Spherical"
            A,b = ForceConSph(x,j)
        elseif j.type == "Weld"
            A,b = ForceConWeldAll(x,j)
        elseif j.type == "Free" || j.type == "Spring"
            A,b = ForceFree(x,j)
        else
            error("Joint type not prescribed.")
        end
    end
    return A,b
end

# For inertial frame and body
function ForceFreeIn(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated due to quaternion constraint
    b2 = j.RB2

    # Revolute Joint has 7 constraints
    A = zeros(T,(1,14));b = T[0.0];

    QcA, Qcb = QuatNormConstraint(x,j)
    A[1,:] = QcA[2,:]; b[1] = Qcb[2]
    return (A,b)

end

function ForceConRevIn(x::AbstractArray{T},j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real#, extF2::extForces, GravityInInertial::Vector{Float64})::Vector{Float64}
    # Constraint Force generated by revolute joint
    b2 = j.RB2

    # Revolute Joint has 6 constraints (wrt inertial frame)
    A = zeros(T,(6,14)); b = zeros(T,6);

    A[1:3,:],b[1:3] = TranslationConstraint(x,j)

    QcA, Qcb = QuatNormConstraint(x,j)
    A[4,:] = QcA[2,:]; b[4] = Qcb[2]

    A[5:end,:] = RevJointConstraint(x,j)[1]; b[5:end] = RevJointConstraint(x,j)[2]

    return (A,b)
end

function ForceConRev2In(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by 2revolute joint

    # Revolute Joint has 5 constraints (wrt inertial frame)
    A = zeros(T,(5,14)); b = zeros(T,5);
    A[1:3,:],b[1:3] = TranslationConstraint(j);

    QcA, Qcb = QuatNormConstraint(x,j)
    A[4,:] = QcA[2,:]; b[4] = Qcb[2]

    A[5,:],b[5]  = RevJointConstraint(j);

    return (A,b)
end

function ForceConWeldAllIn(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by weld joint

    # Weld Joint has 7 constraints
    A = zeros(T,(8,14)); b = zeros(T,8);

    A[1:3,:],b[1:3] = TranslationConstraint(j);

    QcA, Qcb = QuatNormConstraint(x,j);
    A[4,:] = QcA[2,:]; b[4] = Qcb[2];

    A[5:end,:], b[5:end] = WeldJointAllConstraint(j);

    return (A,b)
end

function ForceConSphIn(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by spherical joint

    # Spherical Joint leaves 3 degrees of freedom (remove 4 constraints)
    A = zeros(T,(4,14)); b = zeros(T,4)

    A[1:3,:],b[1:3] = TranslationConstraint(x,j);

    QcA, Qcb = QuatNormConstraint(x,j)
    A[4,:] = QcA[2,:]; b[4] = Qcb[2]

    return (A,b)
end

# For 2 rigid bodies, no inertial frame
function ForceFree(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated due to quaternion constraint
    A,b = QuatNormConstraint(j);
    return (A,b)
end


function ForceConRev(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by revolute joint
    # Revolute Joint imposes 7 constraints between two rigid bodies - 3 translational, 2 rotational, 2 quat norm
    # Resultant dof of the two bodies: 7
    A = zeros(T, (7,14)); b = zeros(T,7);

    A[1:3,:],b[1:3] = TranslationConstraint(x,j)
    A[4:5,:], b[4:5] = QuatNormConstraint(x,j)
    A[6:end,:],b[6:end] = RevJointConstraint(x,j)

    return (A,b)
end

function ForceConRev2(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by revolute2 joint
    # Revolute Joint2 has 6 constraints - 3 translational, 2 quat norm, 1 rotational
    A = zeros(T,(6,14)); b = zeros(T,6);
    A[1:3,:], b[1:3] = TranslationConstraint(j);
    A[4:5,:], b[4:5] = QuatNormConstraint(j);
    A[6,:], b[6] = RevJoint2Constraint(j);

    return (A,b)
end

function ForceConWeldAll(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by weld joint
    # Weld Joint has 8 constraints
    A = zeros(T,(8,14)); b = zeros(T,8);
    A[1:3,:],b[1:3] = TranslationConstraint(x,j)
    A[4:5,:], b[4:5] = QuatNormConstraint(x,j)
    A[6:end,:], b[6:end] = WeldJointAllConstraint(j)

    return (A,b)
end

function ForceConSpr(x::AbstractArray{T}, j::Joint) where T<:Real
    b1_id = j.RB1.bodyID; b2_id = j.RB2.bodyID;
    x1 = x[(b1_id-1)*14+1:b1_id*14]
    x2 = x[(b2_id-1)*14+1:b2_id*14]
    β1 = x1[4:7]; β2 = x2[4:7]
    b1_dcm = quat2dcm(β1); b2_dcm = quat2dcm(β2);

    # To compute Spring force and torque
    posSpr1 = x1[1:3] + transpose(b1_dcm)*j.pos1 # Inertial position of spring connection on first body
    posSpr2 = x2[1:3] + transpose(b2_dcm)*j.pos2 # Inertial position of spring connection on second body
    unitVec = (posSpr2 - posSpr1)/norm(posSpr2-posSpr1,2) # unit Vector in direction of spring (first body to second body)

    # First body
    E1 = genE(β1)
    F1 = j.k*(norm(posSpr1-posSpr2)-j.restLen)*unitVec # Force exerted by spring on the first body (inertial frame)
    F1 += j.c*((x2[8:10] - x1[8:10])); # Damper force (damper coefficient*(velocity difference))
    τ1 = cross(transpose(b1_dcm)*j.pos1,F1) # Torque exerted by spring-damper on the first body (inertial frame)
    Γb1 = [0.0;b1_dcm*τ1] #(body frame)
    Γu1 = 2*transpose(E1)*Γb1 # generalized torque vector

    # Second body
    E2 = genE(β2)
    F2 = -F1 # Force exerted by spring-damper on second body
    τ2 = cross(transpose(b2_dcm)*j.pos2,F2) # Torque exerted by spring on second body
    Γb2 = [0.0;b2_dcm*τ2]
    Γu2 = 2*transpose(E2)*Γb2

    F_b1 = [F1; Γu1]; F_b2 = [F2; Γu2];
    return (F_b1, F_b2)
end

function ForceConSph(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Constraint Force generated by spherical joint
    # Spherical Joint has 5 constraints (3 translational, 2 quat norm)
    A = zeros(T,(5,14)); b = zeros(T,5)
    A[1:3,:], b[1:3] = TranslationConstraint(j);
    A[4:5,:],b[4:5] = QuatNormConstraint(j)

    return (A,b)
end

function TranslationConstraint(x::AbstractArray{T},j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    b1_id = j.RB1.bodyID; b2_id = j.RB2.bodyID;
    x1 = x[(b1_id-1)*14+1:b1_id*14]
    x2 = x[(b2_id-1)*14+1:b2_id*14]
    b1 = j.RB1
    b2 = j.RB2
    rj1 = j.pos1
    rj2 = j.pos2

    A = zeros(T,(3,14)); b = zeros(T,3)
    A[:,1:3] = Matrix{Float64}(I,3,3)
    A[:,8:10] = -Matrix{Float64}(I,3,3)

    b1Tcs = TranslationConstraintSupplement(x1, rj1)
    b2Tcs = TranslationConstraintSupplement(x2, rj2)
    A[:,4:7] = b1Tcs[1]
    A[:,11:14] = -b2Tcs[1]
    b = -b1Tcs[2] + b2Tcs[2]

    return (A,b)

end

function TranslationConstraintSupplement(x::AbstractArray{T},pos::Vector{Float64}) where T <: Real
    # To constrain joint wrt body's cm
    β = x[4:7]
    βdot = x[11:14]

    r(y::AbstractArray{T}) where T<:Real = transpose(quat2dcm(y))*pos
    rJac = z->ForwardDiff.jacobian(r,z)

    rdot = y->ForwardDiff.jacobian(r,y)*βdot
    rβdot = y->ForwardDiff.jacobian(r,β)*y
    rddotRHS = ForwardDiff.jacobian(rdot,β)*βdot
    rddotLHS = rJac(β)
    return (rddotLHS,rddotRHS)
end

function QuatNormConstraint(x::AbstractArray{T},j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    A = zeros(T,(2,14))
    b = zeros(T,2)

    b1 = j.RB1
    b2 = j.RB2

    b1_id = b1.bodyID;
    b2_id = b2.bodyID;
    β1 = x[14*(b1_id-1)+4:14*(b1_id-1)+7]
    β2 = x[14*(b2_id-1)+4:14*(b2_id-1)+7]
    β1dot = x[14*(b1_id-1)+11:14*(b1_id-1)+14]
    β2dot = x[14*(b2_id-1)+11:14*(b2_id-1)+14]

    A[1,4:7] = β1
    A[2,11:14] = β2
    b[1] = -transpose(β1dot)*β1dot
    b[2] = -transpose(β2dot)*β2dot
    return (A,b)
end


function RevJointConstraint(x::AbstractArray{T},j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Attempting relative angular ω between bodies in body 1 frame only has component about the joint axis
    b1 = j.RB1
    b2 = j.RB2
    b1_id = b1.bodyID;
    b2_id = b2.bodyID;
    β1 = x[14*(b1_id-1)+4:14*(b1_id-1)+7]
    β2 = x[14*(b2_id-1)+4:14*(b2_id-1)+7]
    β1dot = x[14*(b1_id-1)+11:14*(b1_id-1)+14]
    β2dot = x[14*(b2_id-1)+11:14*(b2_id-1)+14]
    axis = j.axis
    βaug = [β1;β1dot;β2;β2dot]
    function f(y::AbstractArray{T}) where T <: Real
        y1  = y[1:4]; y1dot = y[5:8]#β1; β1dot
        y2 = y[9:12]; y2dot = y[13:16] #β2;β2dot
        ωb2 = angVel(y2,y2dot)
        ωb1 = angVel(y1,y1dot)
        ωb2inb1 = quat2dcm(y1)*transpose(quat2dcm(y2))*ωb2
        relAngVel = ωb2inb1 - ωb1
        return relAngVel
    end
    fJac = z->ForwardDiff.jacobian(f,z)

    Ax = zeros(T,(3,14)); bx = zeros(T,3)
    Ax[:,4:7] = fJac(βaug)[:,5:8]
    Ax[:,11:14] = fJac(βaug)[:,13:16]
    bx = -fJac(βaug)[:,1:4]*β1dot - fJac(βaug)[:,9:12]*β2dot

    zAxis = findall(iszero,axis)
    A = zeros(T,(2,14)); b = zeros(T,2)
    A[1,:] = Ax[zAxis[1],:]
    A[2,:] = Ax[zAxis[2],:]
    b[1] = bx[zAxis[1]]
    b[2] = bx[zAxis[2]]

    return (A,b)

end

function RevJoint2Constraint(x::AbstractArray{T}, j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    b1 = j.RB1
    b2 = j.RB2
    b1_id = b1.bodyID;
    b2_id = b2.bodyID;
    x1 = x[14*(b1_id-1)+1:14*(b1_id-1)+14];
    x2 = x[14*(b2_id-1)+1:14*(b2_id-1)+14];

    A = zeros(T,(1,14));
    b = T(0.0)
    y = RevJoint2ConstraintSupplement(x1,x2)
    A[:,4:7] = y[1][1]
    A[:,11:14] = y[1][2]
    b = y[2][1]
    return (A,b)
end

function RevJoint2ConstraintSupplement(x1::AbstractArray{T},x2::AbstractArray{T}) where T <: Real
    function f(y::AbstractArray{T}) where T <: Real
        # Computes the Hamilton product of β2 and β1^-1 (i.e. β2*β1^-1))
        # Let Q be the product.
        # (Q(1) + Q(2))^2 + (Q(0) + Q(3))^2 = 1.
        q1 = [y[1];-y[2:4]]
        q2 = y[5:8]
        Q = quaternionProduct(q2,q1)
        # Q = [q2[1]*q1[1] - transpose(q2[2:4])*q1[2:4];
        #      q2[1]*q1[2:4] + q1[1]*q2[2:4] + cross(q2[2:4],q1[2:4])]
        x = (Q[2] + Q[3])^2 + (Q[1] + Q[4])^2
        return x
    end
    fGrad = z->ForwardDiff.gradient(f,z)
    fdotJac = z->ForwardDiff.jacobian(fGrad,z)

    β1 = x1[4:7]
    β2 = x2[4:7]
    β1dot = x1[11:14]
    β2dot = x2[11:14]
    β = [β1;β2]

    y1 = f(β)
    y2 = fGrad(β)
    y3 = fdotJac(β)
    y3r = reshape(y3,(1,8,8))
    x1 = y3[1:4,1:4]*β1dot; x1 = reshape(x1,(1,4))
    x2 = y3[1:4,5:8]*β2dot; x2 = reshape(x2,(1,4))
    x3 = y3[5:end,1:4]*β1dot; x3 = reshape(x3,(1,4))
    x4 = y3[5:end,5:8]*β2dot; x4 = reshape(x4,(1,4))
    # β = rand(8); β1dot = rand(4); β2dot = rand(4);
    LHS = (reshape(-fGrad(β)[1:4],(1,4)),reshape(-fGrad(β)[5:8],(1,4)))
    RHS = (x1+x2)*β1dot + (x3+x4)*β2dot
    return (LHS,RHS)

end

function WeldJointAllConstraint(x::AbstractArray{T},j::Joint)::Tuple{AbstractArray{T}, AbstractArray{T}} where T<:Real
    # Attempting relative angular ω between bodies in body 1 frame is zero
    b1 = j.RB1
    b2 = j.RB2
    b1_id = b1.bodyID;
    b2_id = b2.bodyID;
    x1 = x[14*(b1_id-1)+1:14*(b1_id-1)+14];
    x2 = x[14*(b2_id-1)+1:14*(b2_id-1)+14];
    β1 = x1[4:7]
    β2 = x2[4:7]
    β1dot = x1[11:14]
    β2dot = x2[11:14]
    axis = j.axis
    βaug = [β1;β1dot;β2;β2dot]
    function f(y::AbstractArray{T}) where T <: Real
        y1  = y[1:4]; y1dot = y[5:8]#β1
        y2 = y[9:12]; y2dot = y[13:16] #β2;β2dot
        ωb2 = angVel(y2,y2dot)
        ωb1 = angVel(y1,y1dot)
        ωb2inb1 = quat2dcm(y1)*transpose(quat2dcm(y2))*ωb2
        relAngVel = ωb2inb1 - ωb1
        return relAngVel
    end
    fJac = z->ForwardDiff.jacobian(f,z)

    Ax = zeros(T,(3,14)); bx = zeros(T,3)
    Ax[:,4:7] = fJac(βaug)[:,5:8]
    Ax[:,11:14] = fJac(βaug)[:,13:16]
    bx = -fJac(βaug)[:,1:4]*β1dot - fJac(βaug)[:,9:12]*β2dot

    return (Ax,bx)
end

function genMatM(X::AbstractArray{T},b::RigidBody)::AbstractArray{T} where T<:Real
    # Function to generate mass matrix for each rigid body
    # T <: Real
    b_id = b.bodyID;
    β = X[14*(b_id-1)+4:14*(b_id-1)+7]
    E = genE(β)
    J = b.J
    M = [b.m*Matrix{Float64}(I,3,3)          zeros(3,4)
                         zeros(4,3) 4*transpose(E)*J*E]
    return M
end

# function genExtF(X::AbstractArray{T},b::RigidBody,extF::extForces,GravityInInertial::MArray{Tuple{3},Real,1,3})::AbstractArray{T} where T<:Real
#     # Function to generate augmented external Force vector for unconstrained system
#     # External Forces are always in the body frame
#     b_id = b.bodyID;
#     β = X[14*(b_id-1)+4:14*(b_id-1)+7]
#     βdot = X[14*(b_id-1)+11:14*(b_id-1)+14]
#     dcm = quat2dcm(β)
#     # b.ω = angVel(β,βdot)
#     E = genE(β)
#     Edot = genE(βdot)
#     TotalMoment = zeros(3)
#     for i in 1:size(extF.Forces)[1]
#         TotalMoment = TotalMoment + cross(extF.Positions[i,:],extF.Forces[i,:])
#     end
#     TotalMoment = TotalMoment + sum(extF.Torques,dims=1)[:]
#     Γb = [0.0;TotalMoment] # In the body frame
#     Γu = 2*transpose(E)*Γb
#     F = [transpose(dcm)*(sum(extF.Forces,dims=1)[:]) + b.m*GravityInInertial
#          Γu - 8*transpose(Edot)*b.J*E*βdot - 4*b.J[1,1]*(transpose(βdot)*βdot)*β]
#     return F
#     # F is a 7x1 vector in the inertial frame
# end

function genExtF(X::AbstractArray{T},b::RigidBody,U::AbstractArray{S},GravityInInertial::Vector{Float64})::AbstractArray{Union{T,S}} where {T<:Real, S<:Real}
    # Function to generate augmented external Force vector for unconstrained system
    # External Forces are always in the body frame
    b_id = b.bodyID;
    β = X[14*(b_id-1)+4:14*(b_id-1)+7]
    βdot = X[14*(b_id-1)+11:14*(b_id-1)+14]
    dcm = quat2dcm(β)
    # b.ω = angVel(β,βdot)
    E = genE(β)
    Edot = genE(βdot)
    TotalMoment = U[4:6]

    Γb = [0.0;TotalMoment] # In the body frame
    Γu = 2*transpose(E)*Γb

    F = Vector{Union{T,S}}(undef,7)
    F = [transpose(dcm)*U[1:3] + b.m*GravityInInertial
          Γu - 8*transpose(Edot)*b.J*E*βdot - 4*b.J[1,1]*(transpose(βdot)*βdot)*β]
    return F
    # F is a 7x1 vector in the inertial frame
end

##
function unique_ids(A::AbstractArray{T,2}) where T<:Real
    # Function to remove the redundant rows in constraint matrix A that repeat the quaternion norm constraint for a body connected by multiple joints

    A_unq = unique(A,dims=1); #after removing duplicate rows from A
    inds = Int64[]; # indices of non unique rows
    for i=1:size(A_unq,1)
        for j=1:size(A,1)
            if A_unq[i,:] == A[j,:]
                push!(inds,j)
                break;
            end
        end
    end
    return A_unq, inds
end

##
# trying StaticVector
