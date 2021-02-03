include("../src/simulateDiff.jl")

using ForwardDiff
using FiniteDifferences
m = central_fdm(5,1); # first order differential using fifth order central difference method

function kronSum(A::AbstractArray{T,2},B::AbstractArray{T,2}) where T<:Real
    out = kron(A, Matrix{T}(I,(size(B)))) + kron(Matrix{T}(I,(size(A))),B)
    return out
end

function linearize(x,u,j,GravityInInertial)
    nJ = length(j);

    function MDiff(x::Vector{T}) where T<:Real
        M = assembleM(x,j)
        Ms = real(sqrt(M))
        M_inv = inv(M)
        Ms_inv = real(sqrt(M_inv))

        dM = ForwardDiff.jacobian(z -> assembleM(z,j),x)
        dM_inv = -(kron(permutedims(inv(M)),inv(M)))*ForwardDiff.jacobian(z->assembleM(z,j), x)
        dMs = inv(kronSum(permutedims(Ms), Ms))*dM
        dMs_inv = -kron(permutedims(Ms_inv),Ms_inv)*(dMs);
        return M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv
    end
    # M, dM, M_inv, dM_inv, Ms, dMs, Ms_inv, dMs_inv = MDiff(x);
    # A, b = Ab_VecOfMat(x,j)
    # Fu = assembleF(x,u,j,GravityInInertial)

    function t1Fn(x)
        A, b = Ab_VecOfMat(x,j)
        M = assembleM(x,j)
        M, _, _, _, _, _, M_inv, dM_inv = MDiff(x);
        # M = genMatM(x,j[1].RB2)
        # M1 = M#[4:7,4:7]
        # F = svd(M);
        # M_inv = F.V*diagm(1./(F.S))*(permutedims(F.U))
        # println("err = ", norm(inv(M) - pinv(M)))
        # sleep(1000)
        # M, _, Ms_inv, _, _, dMs_inv = MDiff(x)
        sz_A1 = size(A,1)
        # t1 = permutedims(A)*((A*(M\(permutedims(A))))\I(sz_A1))
        t1 = permutedims(A)*inv(A*inv(M)*permutedims(A))

        dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
        dA_tr = ForwardDiff.jacobian(z->permutedims(Ab_VecOfMat(z,j)[1]),x);

        Y = A*M_inv*permutedims(A); invY = inv(Y);
        dY1 = kron(A*M_inv,1.0I(sz_A1))*dA;
        dY2 = kron(A,A)*dM_inv;
        dY3 = kron(1.0I(sz_A1), A*M_inv)*dA_tr;
        dY = dY1 + dY2 + dY3;
        dY_inv = -kron(permutedims(invY),invY)*dY;

        dt1_1 = kron(permutedims(invY),1.0I(size(A,2)))*dA_tr;
        dt1_2 = kron(1.0I(sz_A1),permutedims(A))*dY_inv;
        dt1 = dt1_1 + dt1_2;
        return t1, dt1
    end

    t1, dt1 = t1Fn(x);
    function GDiff(x::Vector{T}) where T<:Real
        A, b = Ab_VecOfMat(x,j)
        M, _, Ms_inv, _, _, dMs_inv, _, _ = MDiff(x)

        G = A*Ms_inv;
        Gp = permutedims(G)*((G*permutedims(G))\I(size(G,1)))

        t1, dt1 = t1Fn(x)
        dGp1 = kron(permutedims(t1),Matrix{T}(I,size(M)))*dMs_inv

        sz_t2 =size(t1,2)
        t2 = kron(Matrix{T}(I,(sz_t2,sz_t2)), Ms_inv)
        dGp2 = t2*dt1;

        dGp = dGp1 + dGp2
        return Gp, dGp
    end
    # Gp, dGp = GDiff(x);

    function hDiff(x)
        Fu = assembleF(x,u,j,GravityInInertial)
        A, b = Ab_VecOfMat(x,j)
        M = assembleM(x,j)

        a = inv(M)*Fu
        h = b-A*a
        return h
    end
    # h = hDiff(x);

    function FcDiff(x::Vector{T}) where T<:Real
        Fu = assembleF(x,u,j,GravityInInertial)
        A, b = Ab_VecOfMat(x,j)
        M, Ms, _, _, dMs, _, _, _ = MDiff(x)

        a = M\Fu
        h = b-A*a

        Gp, dGp = GDiff(x)

        Fc = Ms*Gp*h

        dFc1 = kron(permutedims(Gp*h),I(size(M,1)))*dMs
        dFc2 = kron(permutedims(h),Ms)*dGp
        dFc3 = Ms*Gp*ForwardDiff.jacobian(z->hDiff(z),x)

        dFc = dFc1 + dFc2 + dFc3;

        dFu = ForwardDiff.jacobian(z->assembleF(z,u,j,GravityInInertial),x)
        return Fc, dFc, Fu, dFu
    end

    # Fc, dFc, Fu, dFu = FcDiff(x)

    function accDiff(x::Vector{T}) where T<:Real
        len_x = length(x)-14

        M = assembleM(x,j)
        Fc, dFc, Fu, dFu = FcDiff(x)

        dacc = zeros(T,(7*nJ,length(x)))
        for k=1:nJ
            invJ = j[k].RB2.invJ;
            ct1 = 7*(k-1); ct2 = 14*k;
            dacc[ct1+1:ct1+3,:] = 1/j[k].RB2.m*(dFc[ct1+1:ct1+3,:] +   dFu[ct1+1:ct1+3,:])


            E = genE(x[ct2+4:ct2+7])
            fdJ_E = ForwardDiff.jacobian(z->genE(z[ct2+4:ct2+7]),x)
            fdJ_trE = ForwardDiff.jacobian(z->transpose(genE(z[ct2+4:ct2+7])),x)

            dacc2_1 = kron(permutedims(invJ*E*((Fc+Fu)[ct1+4:ct1+7])),I(size(E,1)))*fdJ_trE
            dacc2_2 = kron(permutedims((Fc+Fu)[ct1+4:ct1+7]),permutedims(E)*invJ)*fdJ_E
            dacc2_3 = permutedims(E)*invJ*E*(dFc[ct1+4:ct1+7,:] + dFu[ct1+4:ct1+7,:])

            dacc[ct1+4:ct1+7,:] = 1/4*(dacc2_1 + dacc2_2 + dacc2_3)
        end
        return dacc
    end

    # dacc = accDiff(x);

    function fxdotDiff(x::Vector{T}) where T<:Real
        len_x = length(x)
        dFx = zeros(T,(len_x,len_x))
        # dFx = Matrix{T}(undef,(len_x,len_x))

        dacc = accDiff(x);
        # println("dacc = ", dacc)

        for k=1:length(j)
            dFx[14*k+1:14*k+7,14*k+8:14*k+14] = 1.0*I(7)
            dFx[14*k+8:14*k+14,:] = dacc[7*(k-1)+1:7*k,:]
        end
        return dFx
    end

    dFx = fxdotDiff(x)
    dFu = ForwardDiff.jacobian(z -> fxdot(x,z,j,GravityInInertial),u)
 return dFx, dFu
    # return dFc
end
## Test for different robots
x0, u0 = getXU_0(j);
fdJ_A, fdJ_B = linearize(x0,u0,j,g);
finJ_A = FiniteDifferences.jacobian(m, z->fxdot(z,u0,j,g),x0)[1];
finJ_B = FiniteDifferences.jacobian(m, z->fxdot(x0,z,j,g),u0)[1];


println("err_A = ", norm(fdJ_A - finJ_A))
println("err_B = ", norm(fdJ_B - finJ_B))

@btime linearize(x0,u0,j,g);
@btime FiniteDifferences.jacobian(m, z->fxdot(z,u0,j,g),x0);
