function angle2quat(q::Vector{Float64},ang::Vector{Float64},seq::String)
    c = cos.(ang*0.5);
    s = sin.(ang*0.5);

    if seq=="321"
        q[1] = c[1]*c[2]*c[3] + s[1]*s[2]*s[3]; # β0
        q[2] = c[1]*c[2]*s[3] - s[1]*s[2]*c[3]; # β1
        q[3] = c[1]*s[2]*c[3] + s[1]*c[2]*s[3]; # β2
        q[4] = s[1]*c[2]*c[3] - c[1]*s[2]*s[3]; # β3
    end
end

function angle2quat(ang::Vector{Float64},seq::String)::Vector{Float64}
    q = Vector{Float64}(undef,4);
    angle2quat(q,ang,seq);
    return(q);
end

function quat2angle(ang::Vector{Float64},q::Vector{Float64},seq::String)
    if seq == "321" # Have to program other sequences, but we may not use it.
        # roll (x-axis rotation)
        sinr_cosp = +2.0 * (q[1] * q[2] + q[3] * q[4]);
        cosr_cosp = +1.0 - 2.0 * (q[2]^2 + q[3]^2);
        ϕ = atan(sinr_cosp, cosr_cosp);

        # pitch (y-axis rotation)
        sinp = +2.0 * (q[1] * q[3] - q[4] * q[2]);
        if (abs(sinp) >= 1)
            θ = copysign(M_PI / 2, sinp); # use 90 degrees if out of range
        else
            θ = asin(sinp);
        end

        # yaw (z-axis rotation)
        siny_cosp = +2.0 * (q[1] * q[4] + q[2] * q[3]);
        cosy_cosp = +1.0 - 2.0 * (q[3] * q[3] + q[4] * q[4]);
        ψ = atan(siny_cosp, cosy_cosp);

        ang[1] = ψ;
        ang[2] = θ;
        ang[3] = ϕ;
    end
end

function quat2angle(q::Vector{Float64},seq::String)::Vector{Float64}
    ang = Vector{Float64}(undef,3);
    quat2angle(ang,q,seq)
    return(ang);
end

function dcm2quat(q::Vector{Float64},C::Matrix{Float64})
    # Check dimensions, etc.

    # Apply Sheppard's algorithm
    trace = tr(C);

    v0 = 0.25*(1+trace);
    v1 = 0.25*(1+2*C[1,1]-trace);
    v2 = 0.25*(1+2*C[2,2]-trace);
    v3 = 0.25*(1+2*C[3,3]-trace);

    maxV = max(v0,v1,v2,v3); # Take the maximum.
    if (v0 >= maxV)
        q[1] = sqrt(maxV);
        q[2] = (C[2,3]-C[3,2])/4.0/q[1];
        q[3] = (C[3,1]-C[1,3])/4.0/q[1];
        q[4] = (C[1,2]-C[2,1])/4.0/q[1];
    elseif (v1 >= maxV)
        q[2] = sqrt(maxV);
        q[1] = (C[2,3]-C[3,2])/4.0/q[2];
        q[3] = (C[1,2]+C[2,1])/4.0/q[2];
        q[4] = (C[3,1]+C[1,3])/4.0/q[2];
    elseif (v2 >= maxV)
        q[3] = sqrt(maxV);
        q[1] = (C[3,1]-C[1,3])/4.0/q[3];
        q[2] = (C[1,2]+C[2,1])/4.0/q[3];
        q[4] = (C[2,3]+C[3,2])/4.0/q[3];
    else
        q[4] = sqrt(maxV);
        q[1] = (C[1,2]-C[2,1])/4.0/q[4];
        q[2] = (C[3,1]+C[1,3])/4.0/q[4];
        q[3] = (C[2,3]+C[3,2])/4.0/q[4];
    end
end

function dcm2quat(C::Matrix{T}) where T <: Real
    q = Vector{Real}(undef,4)
    # Apply Sheppard's algorithm
    trace = tr(C);

    v0 = 0.25*(1+trace);
    v1 = 0.25*(1+2*C[1,1]-trace);
    v2 = 0.25*(1+2*C[2,2]-trace);
    v3 = 0.25*(1+2*C[3,3]-trace);

    maxV = max(v0,v1,v2,v3); # Take the maximum.
    if (v0 >= maxV)
        q[1] = sqrt(maxV);
        q[2] = (C[2,3]-C[3,2])/4.0/q[1];
        q[3] = (C[3,1]-C[1,3])/4.0/q[1];
        q[4] = (C[1,2]-C[2,1])/4.0/q[1];
    elseif (v1 >= maxV)
        q[2] = sqrt(maxV);
        q[1] = (C[2,3]-C[3,2])/4.0/q[2];
        q[3] = (C[1,2]+C[2,1])/4.0/q[2];
        q[4] = (C[3,1]+C[1,3])/4.0/q[2];
    elseif (v2 >= maxV)
        q[3] = sqrt(maxV);
        q[1] = (C[3,1]-C[1,3])/4.0/q[3];
        q[2] = (C[1,2]+C[2,1])/4.0/q[3];
        q[4] = (C[2,3]+C[3,2])/4.0/q[3];
    else
        q[4] = sqrt(maxV);
        q[1] = (C[1,2]-C[2,1])/4.0/q[4];
        q[2] = (C[3,1]+C[1,3])/4.0/q[4];
        q[3] = (C[2,3]+C[3,2])/4.0/q[4];
    end
    return q
end

function quat2dcm(C::Matrix{T},q::Vector{T}) where T <: Real
    β0 = q[1];
    β1 = q[2];
    β2 = q[3];
    β3 = q[4];

    C =  [(β0^2+β1^2-β2^2-β3^2) 2*(β1*β2+β0*β3)           2*(β1*β3-β0*β2);
    2*(β1*β2-β0*β3)         (β0^2-β1^2+β2^2-β3^2)   2*(β2*β3+β0*β1);
    2*(β1*β3+β0*β2)       2*(β2*β3-β0*β1)         (β0^2-β1^2-β2^2+β3^2)];
end

function quat2dcm(q::Vector{T}) where T <: Real
# ForwardDiff has problems with inplace function definitions
    β0 = q[1];
    β1 = q[2];
    β2 = q[3];
    β3 = q[4];

    C =  [(β0^2+β1^2-β2^2-β3^2) 2*(β1*β2+β0*β3)           2*(β1*β3-β0*β2);
    2*(β1*β2-β0*β3)         (β0^2-β1^2+β2^2-β3^2)   2*(β2*β3+β0*β1);
    2*(β1*β3+β0*β2)       2*(β2*β3-β0*β1)         (β0^2-β1^2-β2^2+β3^2)];

    # C = Matrix{Float64}(undef,3,3);
    # quat2dcm(C,q);
    return(C);
end

function ang2dcm(dcm::Matrix{Float64},ang::Vector{Float64},orientation::String)

    θ1 = ang[1]; # phi
    θ2 = ang[2]; # theta
    θ3 = ang[3]; # psi

    c1 = cos(θ1); s1 = sin(θ1);
    c2 = cos(θ2); s2 = sin(θ2);
    c3 = cos(θ3); s3 = sin(θ3);

    M1 = [1 0 0;
          0 c1 s1;
          0 -s1 c1];

    M2 = [c2 0 -s2;
           0 1 0;
          s2 0 c2];

    M3 = [c3 s3 0;
         -s3 c3 0;
           0  1 1];

    if orientation == "321"
        C = M1*M2*M3; # Inertial to body
    else
        error("Uknown orientation specified")
    end
end

function skewX(x::Vector{T}) where T<: Real
    X = [    0 -x[3]  x[2]
          x[3]     0 -x[1]
         -x[2]  x[1]    0]
    return X
end

function angVel(β::Vector{T},βdot::Vector{T}) where T<: Real
    β0 = β[1];
    β1 = β[2];
    β2 = β[3];
    β3 = β[4];

    E1 = [-β1  β0  β3 -β2
          -β2 -β3  β0  β1
          -β3  β2 -β1  β0]
    ω = 2*E1*βdot
    return ω
end

function genE(β::Vector{Float64})
    X = [                        transpose(β)
         -β[2:4] β[1]*Matrix{Float64}(I,3,3) - skewX(β[2:4])]
    return X

end

# function quatProdRev(β1::Vector{Float64},β2::Vector{Float64})
#     # Returns the Hamilton product β2*β1^-1 for the joint transformation between bodies 1 and 2
#     β1inv = [β1[1];-β1[2:4]]
#     β2[2:4] = quat2dcm(β1inv)*β2[2:4]
#     x = [β2[1]*β1[1] - transpose(β2[2:4])*(-β1[2:4]);
#          β2[1]*(-β1[2:4]) + β1[1]*β2[2:4] + cross(β2[2:4],-β1[2:4])]
#     return x
# end
#
# function axixRev(β1::Vector{Float64},β2::Vector{Float64})
#     y = quatProdRev(β1,β2)
#     axisInertial= y[2:4]/norm(y[2:4],2)
#     return axisBody = quat2dcm(β1)*axisInertial
# end

function quaternionProduct(β1::Vector{T},β2::Vector{T}) where T <: Real
    x = [β1[1]*β2[1] - dot(β1[2:4],β2[2:4]);
         β1[2]*β2[2:4] + β2[1]*β1[2:4] + cross(β1[2:4],β2[2:4])]
    return x
end
