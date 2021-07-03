J = j[1].RB2.J;
β2 = j[1].RB2.x[4:7]
# β2 = []
function M_β(β::Vector{T}) where T<:Real
    E = genE(β)
    M = 4*transpose(E)*J*E
    if typeof(β[1]) == Float64
        Ms = real(sqrt(M))
        # C = cholesky(M)
        # println("L = ", norm(Ms - C.L))
        # println(norm(Ms - C.U))
        # Ms = (sqrt.(diag(M))).*Matrix{T}(I,size(M))
    else # ForwardDiff.Dual variable
        # C = cholesky(M)
        # Ms = C.U

        # Ms = (sqrt.(diag(M))).*Matrix{T}(I,size(M))

        # diag_Mfinal = diag(Mfinal);
        # vals_sq = zeros(T,size(Mfinal,1))
        # for i=1:size(Mfinal,1)
        #     vals_sq[i] = sqrt(diag_Mfinal[i])
        # end
        # Ms = diagm(vals_sq);

        # Ms = (Mfinal^(0.5))

        vals, vecs = eigen(M);
        vals_sq = zeros(T,size(vals))
        for i=1:length(vals)
            vals_sq[i] = sqrt(vals[i])
        end
        Ms = vecs*diagm(vals_sq)*(vecs\I(size(vecs,1)))
    end
    return M, Ms
end
# M_β(β2)
m = central_fdm(10,1)
temp1 = FiniteDifferences.jacobian(m, z->M_β(z)[2], β2)[1]
temp2 = ForwardDiff.jacobian(z->M_β(z)[1], β2)
temp1_nz = findall(x -> abs(x)>1e-6, temp1)
temp2_nz = findall(x -> abs(x)>1e-6, temp2)
norm(temp1 - temp2)

# function kronSum(A::Array{T,2},B::Array{T,2}) where T<:Real
#     out = kron(A, Matrix{T}(I,(size(B)))) + kron(Matrix{T}(I,(size(A))),B)
#     return out
# end

M_check, Mhalf_check = M_β(β2)
dM_half = inv(kronSum(permutedims(Mhalf_check), Mhalf_check))
norm(dM_half*temp2 - temp1)

## trying svd for square root
# clearconsole()
# using  FiniteDifferences
# using LinearAlgebra
# using GenericLinearAlgebra
# using ForwardDiff
# using Revise
# function sqrtM(x::Array{T}) where T<:Real
#     # M = assembleM(x,j)
#     M = genMatM(x,j[1].RB2)
#     M1 = M[4:7,4:7]#[5:7,5:7]
#     # E = genE(q)
#     # M = 4*transpose(E)*J*E
#     if typeof(M[1]) == Float64
#         Ms = real(sqrt(M1))
#         # println("Ms = ", Ms)
#         # F = GenericLinearAlgebra.svd!(M)
#         # Ms2 = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
#         # println("Mserr = ", norm(Ms - Ms2))
#     else
#         F = GenericLinearAlgebra.svd!(M1)
#         Ms = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
#     end
#     return M, Ms
# end
# function sqrtMq(q::Array{T}) where T<:Real
#     E = genE(q)
#     M = 4*transpose(E)*j[1].RB2.J*E
#     if typeof(M[1]) == Float64
#         out = real(sqrt(M))
#         # F = GenericLinearAlgebra.svd!(M)
#         # Ms = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
#         # println("Mserr = ", norm(out - Ms))
#     else
#         F = GenericLinearAlgebra.svd!(M)
#         out = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
#     end
#     return M, out
# end
# # x0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.24999999999999997, 0.0, -0.43301270189221935, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.294735414289422];
# # q2 = [sqrt(2)/2; sqrt(2)/2;zeros(2)]
# x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; [1;zeros(3)]; zeros(7)];
# # x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; randomQuat(); zeros(7)];
# updateRigidBody!(j[1].RB2,x0)
# # x0 = randomQuat();
# _,tempMs = sqrtM(x0)
# # temp = rand(10)
# m = central_fdm(10,1)
# finJ = FiniteDifferences.jacobian(m, z->sqrtM(z)[2], x0)[1]
# fdJ = ForwardDiff.jacobian(z->sqrtM(z)[2],x0)
# finJ_nz = findall(x -> abs(x)>1e-6, finJ)
# fdJ_nz = findall( x-> abs(x)>0, fdJ)
# println()
# println(norm(fdJ - finJ))
#
# # q2 = randomQuat();
# q2 = (x0[18:21])
# sqrtMq(q2)
# finJ2 = FiniteDifferences.jacobian(m, z->sqrtMq(z)[2], q2)[1]
# fdJ2 = ForwardDiff.jacobian(z->sqrtMq(z)[2],q2)
# println(norm(fdJ2 - finJ2))
##
# J = j[1].RB2.J;
# function pinvM(q::Array{T}) where T<:Real
#     E = genE(q)
#     M = 4*transpose(E)*J*E
#     Ms = sqrt(diagm(M))
#     F = GenericLinearAlgebra.svd!(diagm(M))
#     out = pinv(diagm(M))
#     return out
# end
# ForwardDiff.jacobian(pinvM, temp)
## symbolic
# using SymPy
#
# u0, u1, u2, u3 = symbols("u0 u1 u2 u3")
#
# E = [u0 u1 u2 u3;
#     -u1 u0 u3 -u2;
#     -u2 -u3 u0 u1;
#     -u3 u2 -u1 u0];
#
# J = j[1].RB2.J;
# M = 4*permutedims(E)*J*E
# F = GenericLinearAlgebra.svd!(M)
# # Ms = (M)^(2)
# diffM = [M.diff((u0,1))[:] M.diff((u1,1))[:] M.diff((u2,1))[:] M.diff((u3,1))[:]]
#
# diffM_eval = zeros(size(diffM))
# for i=1:length(diffM)
#     diffM_eval[i] = diffM[i](1.0, 0.0, 0.0, 0.0)
# end
# diffM_eval
# M[1,1](1.0, 0.0, 0.0, 0.0)
##
# function invM(x)
#     M = diagm(x)
#     return inv(M)
# end
# r=rand(3)
# finJ_M = FiniteDifferences.jacobian(m,z->invM(z), r)[1]
# fdJ_M = ForwardDiff.jacobian(z->invM(z), r)
# println("err= ", norm(finJ_M - fdJ_M))
## Mdiff function
clearconsole()
function kronSum(A::Array{T,2},B::Array{T,2}) where T<:Real
    out = kron(A, Matrix{T}(I,(size(B)))) + kron(Matrix{T}(I,(size(A))),B)
    return out
end
# x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; [1;zeros(3)]; zeros(7)];
# x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; randomQuat(); zeros(7)];

function Mdiff(x::Vector{T}) where T<:Real
    M = assembleM(x,j)
    Ms = real(sqrt(M))
    M_inv = inv(M)
    Ms_inv = real(sqrt(M_inv))

    dM = ForwardDiff.jacobian(z -> assembleM(z,j),x)
    dM_inv = -(kron(permutedims(inv(M)),inv(M)))*ForwardDiff.jacobian(z->assembleM(z,j), x)
    # dM_inv = ForwardDiff.jacobian(z -> inv(assembleM(z,j)),x)
    dMs = inv(kronSum(permutedims(Ms), Ms))*dM
    dMs_inv = -kron(permutedims(Ms_inv),Ms_inv)*(dMs);
    return M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv
end
m = central_fdm(5,1);
# finJ_dMs = FiniteDifferences.jacobian(m, z->Mdiff(z)[2], x0)[1]
# finJ_dMs_inv = FiniteDifferences.jacobian(m, z->Mdiff(z)[3], x0)[1]
# println()
# println(norm(finJ_dMs - Mdiff(x0)[5]))
# println(norm(finJ_dMs_inv - Mdiff(x0)[6]))
#
# _,_,Ms_inv,_,dMs, fdJ_dMs_inv = Mdiff(x0);
# # fdJ_dMs_inv = -kron(permutedims(Ms_inv),Ms_inv)*(dMs);
# println(norm(finJ_dMs_inv - fdJ_dMs_inv))
##
clearconsole()
function t1Fn(x)
    A, b = Ab_VecOfMat(x,j)
    M = assembleM(x,j)
    M, _, _, _, _, _, M_inv, dM_inv = Mdiff(x);
    # M = genMatM(x,j[1].RB2)
    # M1 = M#[4:7,4:7]
    # F = svd(M);
    # M_inv = F.V*diagm(1./(F.S))*(permutedims(F.U))
    # println("err = ", norm(inv(M) - pinv(M)))
    # sleep(1000)
    # M, _, Ms_inv, _, _, dMs_inv = Mdiff(x)
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
function Gdiff(x::Vector{T}) where T<:Real
    A, b = Ab_VecOfMat(x,j)
    M, _, Ms_inv, _, _, dMs_inv, _, _ = Mdiff(x)

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
# Gdiff(x0)
finJ_G = FiniteDifferences.jacobian(m,z->Gdiff(z)[1],x0)[1]
fdJ_G = Gdiff(x0)[2]#ForwardDiff.jacobian(z->Gdiff(z)[2],x0)
println("\nerr_Gp", norm(finJ_G - fdJ_G))
#
# finJ_t1 = FiniteDifferences.jacobian(m,z->t1Fn(z)[1], x0)[1]
# fdJ_t1 = t1Fn(x0)[2]
# println("\nerr_t1= ", norm(finJ_t1 - fdJ_t1))

##
x0, u = getXU_0(j);

function hdiff(x)
    Fu = assembleF(x,u,j,GravityInInertial)
    A, b = Ab_VecOfMat(x,j)
    M = assembleM(x,j)

    a = inv(M)*Fu
    h = b-A*a
    return h
end
function Fcdiff(x::Vector{T}) where T<:Real
    Fu = assembleF(x,u,j,GravityInInertial)
    A, b = Ab_VecOfMat(x,j)
    M, Ms, _, _, dMs, _, _, _ = Mdiff(x)

    a = M\Fu
    h = b-A*a

    Gp, dGp = Gdiff(x)

    Fc = Ms*Gp*h

    dFc1 = kron(permutedims(Gp*h),I(size(M,1)))*dMs
    dFc2 = kron(permutedims(h),Ms)*dGp
    dFc3 = Ms*Gp*ForwardDiff.jacobian(z->hdiff(z),x)

    dFc = dFc1 + dFc2 + dFc3;

    dFu = ForwardDiff.jacobian(z->assembleF(z,u,j,GravityInInertial),x)
    return Fc, dFc, Fu, dFu
end
finJ_Fc = FiniteDifferences.jacobian(m,z->Fcdiff(z)[1],x0)[1]
fdJ_Fc = Fcdiff(x0)[2]
println("\nerr_Fc = ", norm(finJ_Fc - fdJ_Fc))

finJ_Fu = FiniteDifferences.jacobian(m,z->Fcdiff(z)[3],x0)[1]
fdJ_Fu = Fcdiff(x0)[4]
println("\nerr_Fu = ", norm(finJ_Fu - fdJ_Fu))

## fxdotdiff
invJ = j[1].RB2.invJ;
nJ = length(j); # also equal to number of bodies (excluding inertial frame)
function accdiff(x::Vector{T}) where T<:Real
    len_x = length(x)-14
    # dacc = zeros(T,(len_x,len_x))

    M = assembleM(x,j)
    Fc, dFc, Fu, dFu = Fcdiff(x)

    # dM_inv = ForwardDiff.jacobian(z ->inv(assembleM(z,j)), x)

    # dacc1 = kron(permutedims(Fc+Fu),1.0I(size(M,1)))*dM_inv
    # dacc2 = inv(M)*(dFc + dFu)
    #
    # dacc = dacc1 + dacc2
    dacc = zeros(T,(7*nJ,length(x)))
    for k=1:nJ
        ct1 = 7*(k-1); ct2 = 14*k;
        dacc[ct1+1:ct1+3,:] = 1/j[k].RB2.m*(dFc[ct1+1:ct1+3,:] +   dFu[ct1+1:ct1+3,:])


        E = genE(x[ct2+4:ct2+7])
        fdJ_E = ForwardDiff.jacobian(z->genE(z[ct2+4:ct2+7]),x)
        fdJ_trE = ForwardDiff.jacobian(z->transpose(genE(z[ct2+4:ct2+7])),x)

        dacc2_1 = kron(permutedims(invJ*E*((Fc+Fu)[ct1+4:ct1+7])),I(size(E,1)))*fdJ_trE
        dacc2_2 = kron(permutedims((Fc+Fu)[ct1+4:ct1+7]),permutedims(E)*invJ)*fdJ_E
        dacc2_3 = permutedims(E)*invJ*E*(dFc[ct1+4:ct1+7,:] + dFu[ct1+4:ct1+7,:])

        # println("\ndacc2_1 = ", dacc2_1)
        # println("\ndacc2_2 = ", dacc2_2)
        # println("\ndacc2_3 = ", dacc2_3)
        # sleep(1000);
        dacc[ct1+4:ct1+7,:] = 1/4*(dacc2_1 + dacc2_2 + dacc2_3)
    end
    return dacc
end
# accdiff(x0)
##
function fxdotdiff(x::Vector{T}) where T<:Real
    len_x = length(x)
    dFx = zeros(T,(len_x,len_x))
    # dFx = Matrix{T}(undef,(len_x,len_x))

    dacc = accdiff(x);
    # println("dacc = ", dacc)

    for k=1:length(j)
        dFx[14*k+1:14*k+7,14*k+8:14*k+14] = 1.0*I(7)
        dFx[14*k+8:14*k+14,:] = dacc[7*(k-1)+1:7*k,:]
    end
    return dFx
end

finJ_fx = FiniteDifferences.jacobian(m,z->fxdot(z,u,j,GravityInInertial),x0)[1];

fdJ_fx = fxdotdiff(x0);
# println("\nerr_fx = ", norm(finJ_fx[25:end,:] - fdJ_fx[25:end,:]))
println("\nerr_fx = ", norm(finJ_fx- fdJ_fx))
## M for β = [1.0;zeros(3)]
q = deepcopy(x0); q[18:21] = [1.0;zeros(3)]
function Mq(x::AbstractVector{T}) where T<:Real
    M = assembleM(x,j)
    M1 = M[4:7,4:7]
    out = (M1)#\I(size(M1,1))
    return out
end
finJ_M = FiniteDifferences.jacobian(m,z->inv(Mq(z)), q)[1]
M1 = Mq(q);
fdJ_M = -(kron(permutedims(inv(M1)),inv(M1)))*ForwardDiff.jacobian(Mq, q)
println("\nerr_invM = ", norm(finJ_M - fdJ_M))
## symbolic trying in login node (hprc)
# using SymPy
# using LinearAlgebra
# u0, u1, u2, u3 = symbols("u0 u1 u2 u3")
# uSym = [u0;u1;u2;u3];
# E = [u0 u1 u2 u3;
#     -u1 u0 u3 -u2;
#     -u2 -u3 u0 u1;
#     -u3 u2 -u1 u0];
#
# J = j[1].RB2.J;
# M = 4*permutedims(E)*J*E
# invM = (M)^(-1)s
## forwarddiff benchmarks
# function tempFn(z)
#     out = z.^2;
#     return out
# end
# y = rand(100);
# fdJ_fn = x->ForwardDiff.jacobian(tempFn,x)
# @btime fdJ_fn(y);
# @btime ForwardDiff.jacobian(tempFn,y);

jb1_nz = findall(x -> abs(x) > 1e-6, jb1)
jb2_nz = findall(x -> abs(x) > 0, jb2)
