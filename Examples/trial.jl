# # J = j[1].RB2.J;
# # β2 = j[1].RB2.x[4:7]
# # # β2 = []
# # function M_β(β::Vector{T}) where T<:Real
# #     E = genE(β)
# #     M = 4*transpose(E)*J*E
# #     if typeof(β[1]) == Float64
# #         Ms = real(sqrt(M))
# #         # C = cholesky(M)
# #         # println("L = ", norm(Ms - C.L))
# #         # println(norm(Ms - C.U))
# #         # Ms = (sqrt.(diag(M))).*Matrix{T}(I,size(M))
# #     else # ForwardDiff.Dual variable
# #         # C = cholesky(M)
# #         # Ms = C.U
# #
# #         # Ms = (sqrt.(diag(M))).*Matrix{T}(I,size(M))
# #
# #         # diag_Mfinal = diag(Mfinal);
# #         # vals_sq = zeros(T,size(Mfinal,1))
# #         # for i=1:size(Mfinal,1)
# #         #     vals_sq[i] = sqrt(diag_Mfinal[i])
# #         # end
# #         # Ms = diagm(vals_sq);
# #
# #         # Ms = (Mfinal^(0.5))
# #
# #         vals, vecs = eigen(M);
# #         vals_sq = zeros(T,size(vals))
# #         for i=1:length(vals)
# #             vals_sq[i] = sqrt(vals[i])
# #         end
# #         Ms = vecs*diagm(vals_sq)*(vecs\I(size(vecs,1)))
# #     end
# #     return M, Ms
# # end
# # # M_β(β2)
# # m = central_fdm(10,1)
# # temp1 = FiniteDifferences.jacobian(m, z->M_β(z)[2], β2)[1]
# # temp2 = ForwardDiff.jacobian(z->M_β(z)[1], β2)
# # temp1_nz = findall(x -> abs(x)>1e-6, temp1)
# # temp2_nz = findall(x -> abs(x)>1e-6, temp2)
# # norm(temp1 - temp2)
# #
# # # function kronSum(A::Array{T,2},B::Array{T,2}) where T<:Real
# # #     out = kron(A, Matrix{T}(I,(size(B)))) + kron(Matrix{T}(I,(size(A))),B)
# # #     return out
# # # end
# #
# # M_check, Mhalf_check = M_β(β2)
# # dM_half = inv(kronSum(permutedims(Mhalf_check), Mhalf_check))
# # norm(dM_half*temp2 - temp1)
# ## tempFn
# # function tempFn(p,q)
# #     f(p1,q1) = sum(p1.^2) + sum(q1.^2);
# #     df_p = ForwardDiff.gradient(z->f(z,q),p);
# #     return df_p
# # end
# # r1 = rand(3); r2 = rand(3);
# # df_p = tempFn(r1,r2)
# # df_p_u = ForwardDiff.jacobian(z->tempFn(r1,z),r2)
# ## trying svd for square root
# # clearconsole()
# # using  FiniteDifferences
# # using LinearAlgebra
# # using GenericLinearAlgebra
# # using ForwardDiff
# # using Revise
# # function sqrtM(x::Array{T}) where T<:Real
# #     # M = assembleM(x,j)
# #     M = genMatM(x,j[1].RB2)
# #     M1 = M[4:7,4:7]#[5:7,5:7]
# #     # E = genE(q)
# #     # M = 4*transpose(E)*J*E
# #     if typeof(M[1]) == Float64
# #         Ms = real(sqrt(M1))
# #         # println("Ms = ", Ms)
# #         # F = GenericLinearAlgebra.svd!(M)
# #         # Ms2 = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
# #         # println("Mserr = ", norm(Ms - Ms2))
# #     else
# #         F = GenericLinearAlgebra.svd!(M1)
# #         Ms = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
# #     end
# #     return M, Ms
# # end
# # function sqrtMq(q::Array{T}) where T<:Real
# #     E = genE(q)
# #     M = 4*transpose(E)*j[1].RB2.J*E
# #     if typeof(M[1]) == Float64
# #         out = real(sqrt(M))
# #         # F = GenericLinearAlgebra.svd!(M)
# #         # Ms = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
# #         # println("Mserr = ", norm(out - Ms))
# #     else
# #         F = GenericLinearAlgebra.svd!(M)
# #         out = permutedims(F.Vt)*diagm(sqrt.(F.S))*F.Vt
# #     end
# #     return M, out
# # end
# # # x0 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.24999999999999997, 0.0, -0.43301270189221935, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.294735414289422];
# # # q2 = [sqrt(2)/2; sqrt(2)/2;zeros(2)]
# # x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; [1;zeros(3)]; zeros(7)];
# # # x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; randomQuat(); zeros(7)];
# # updateRigidBody!(j[1].RB2,x0)
# # # x0 = randomQuat();
# # _,tempMs = sqrtM(x0)
# # # temp = rand(10)
# # m = central_fdm(10,1)
# # finJ = FiniteDifferences.jacobian(m, z->sqrtM(z)[2], x0)[1]
# # fdJ = ForwardDiff.jacobian(z->sqrtM(z)[2],x0)
# # finJ_nz = findall(x -> abs(x)>1e-6, finJ)
# # fdJ_nz = findall( x-> abs(x)>0, fdJ)
# # println()
# # println(norm(fdJ - finJ))
# #
# # # q2 = randomQuat();
# # q2 = (x0[18:21])
# # sqrtMq(q2)
# # finJ2 = FiniteDifferences.jacobian(m, z->sqrtMq(z)[2], q2)[1]
# # fdJ2 = ForwardDiff.jacobian(z->sqrtMq(z)[2],q2)
# # println(norm(fdJ2 - finJ2))
# ##
# # J = j[1].RB2.J;
# # function pinvM(q::Array{T}) where T<:Real
# #     E = genE(q)
# #     M = 4*transpose(E)*J*E
# #     Ms = sqrt(diagm(M))
# #     F = GenericLinearAlgebra.svd!(diagm(M))
# #     out = pinv(diagm(M))
# #     return out
# # end
# # ForwardDiff.jacobian(pinvM, temp)
# ## symbolic
# # using SymPy
# #
# # u0, u1, u2, u3 = symbols("u0 u1 u2 u3")
# #
# # E = [u0 u1 u2 u3;
# #     -u1 u0 u3 -u2;
# #     -u2 -u3 u0 u1;
# #     -u3 u2 -u1 u0];
# #
# # J = j[1].RB2.J;
# # M = 4*permutedims(E)*J*E
# # F = GenericLinearAlgebra.svd!(M)
# # # Ms = (M)^(2)
# # diffM = [M.diff((u0,1))[:] M.diff((u1,1))[:] M.diff((u2,1))[:] M.diff((u3,1))[:]]
# #
# # diffM_eval = zeros(size(diffM))
# # for i=1:length(diffM)
# #     diffM_eval[i] = diffM[i](1.0, 0.0, 0.0, 0.0)
# # end
# # diffM_eval
# # M[1,1](1.0, 0.0, 0.0, 0.0)
# ##
# # function invM(x)
# #     M = diagm(x)
# #     return inv(M)
# # end
# # r=rand(3)
# # finJ_M = FiniteDifferences.jacobian(m,z->invM(z), r)[1]
# # fdJ_M = ForwardDiff.jacobian(z->invM(z), r)
# # println("err= ", norm(finJ_M - fdJ_M))
# ## Mdiff function
# clearconsole()
# function kronSum(A::Array{T,2},B::Array{T,2}) where T<:Real
#     out = kron(A, Matrix{T}(I,(size(B)))) + kron(Matrix{T}(I,(size(A))),B)
#     return out
# end
# # x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; [1;zeros(3)]; zeros(7)];
# # x0 = [0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.24999999999999997; 0.0; -0.43301270189221935; randomQuat(); zeros(7)];
#
# function Mdiff(x::Vector{T}) where T<:Real
#     M = assembleM(x,j)
#     Ms = real(sqrt(M))
#     M_inv = inv(M)
#     Ms_inv = real(sqrt(M_inv))
#
#     dM = ForwardDiff.jacobian(z -> assembleM(z,j),x)
#     dM_inv = -(kron(permutedims(M_inv),M_inv))*dM;
#     dMs = inv(kronSum(permutedims(Ms), Ms))*dM
#     dMs_inv = -kron(permutedims(Ms_inv),Ms_inv)*(dMs);
#     return M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv
# end
# m = central_fdm(5,1);
# # finJ_M_inv = FiniteDifferences.jacobian(m, z->inv(assembleM(z,j)), x0)[1];
# # finJ_M = FiniteDifferences.jacobian(m, z->(assembleM(z,j)), x0)[1];
# # norm(dM_inv - finJ_M_inv)
# # norm(dM - finJ_M)
# # finJ_dMs = FiniteDifferences.jacobian(m, z->Mdiff(z)[2], x0)[1]
# # finJ_dMs_inv = FiniteDifferences.jacobian(m, z->Mdiff(z)[3], x0)[1]
# # println()
# # println(norm(finJ_dMs - Mdiff(x0)[5]))
# # println(norm(finJ_dMs_inv - Mdiff(x0)[6]))
# #
# # _,_,Ms_inv,_,dMs, fdJ_dMs_inv = Mdiff(x0);
# # # fdJ_dMs_inv = -kron(permutedims(Ms_inv),Ms_inv)*(dMs);
# # println(norm(finJ_dMs_inv - fdJ_dMs_inv))
# ##
# clearconsole()
# function t1Fn(x)
#     A, b = Ab_VecOfMat(x,j)
#     M = assembleM(x,j)
#     M, _, _, _, _, _, M_inv, dM_inv = Mdiff(x);
#     # M = genMatM(x,j[1].RB2)
#     # M1 = M#[4:7,4:7]
#     # F = svd(M);
#     # M_inv = F.V*diagm(1./(F.S))*(permutedims(F.U))
#     # println("err = ", norm(inv(M) - pinv(M)))
#     # sleep(1000)
#     # M, _, Ms_inv, _, _, dMs_inv = Mdiff(x)
#     sz_A1 = size(A,1)
#     # t1 = permutedims(A)*((A*(M\(permutedims(A))))\I(sz_A1))
#     t1 = permutedims(A)*inv(A*inv(M)*permutedims(A))
#
#     dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
#     dA_tr = ForwardDiff.jacobian(z->permutedims(Ab_VecOfMat(z,j)[1]),x);
#
#     Y = A*M_inv*permutedims(A); invY = inv(Y);
#     dY1 = kron(A*M_inv,1.0I(sz_A1))*dA;
#     dY2 = kron(A,A)*dM_inv;
#     dY3 = kron(1.0I(sz_A1), A*M_inv)*dA_tr;
#     dY = dY1 + dY2 + dY3;
#     dY_inv = -kron(permutedims(invY),invY)*dY;
#
#     dt1_1 = kron(permutedims(invY),1.0I(size(A,2)))*dA_tr;
#     dt1_2 = kron(1.0I(sz_A1),permutedims(A))*dY_inv;
#     dt1 = dt1_1 + dt1_2;
#     return t1, dt1
# end
# function Gdiff(x::Vector{T}) where T<:Real
#     A, b = Ab_VecOfMat(x,j)
#     M, _, Ms_inv, _, _, dMs_inv, _, _ = Mdiff(x)
#
#     G = A*Ms_inv;
#     Gp = permutedims(G)*((G*permutedims(G))\I(size(G,1)))
#
#     t1, dt1 = t1Fn(x)
#     dGp1 = kron(permutedims(t1),Matrix{T}(I,size(M)))*dMs_inv
#
#     sz_t2 =size(t1,2)
#     t2 = kron(Matrix{T}(I,(sz_t2,sz_t2)), Ms_inv)
#     dGp2 = t2*dt1;
#
#     dGp = dGp1 + dGp2
#     return Gp, dGp
# end
# # Gdiff(x0)
# # finJ_G = FiniteDifferences.jacobian(m,z->Gdiff(z)[1],x0)[1]
# # fdJ_G = Gdiff(x0)[2]#ForwardDiff.jacobian(z->Gdiff(z)[2],x0)
# # println("\nerr_Gp", norm(finJ_G - fdJ_G))
# #
# # finJ_t1 = FiniteDifferences.jacobian(m,z->t1Fn(z)[1], x0)[1]
# # fdJ_t1 = t1Fn(x0)[2]
# # println("\nerr_t1= ", norm(finJ_t1 - fdJ_t1))
#
# ##
# x0, u = getXU_0(j);
# M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x0);
#
# function hdiff(x)
#     Fu = assembleF(x,u,j,GravityInInertial)
#     A, b = Ab_VecOfMat(x,j)
#     M = assembleM(x,j)
#
#     a = inv(M)*Fu
#     h = b-A*a
#     return h
# end
# function Fcdiff(x::Vector{T}) where T<:Real
#     Fu = assembleF(x,u,j,GravityInInertial)
#     A, b = Ab_VecOfMat(x,j)
#     M, Ms, _, _, dMs, _, _, _ = Mdiff(x)
#
#     a = M\Fu
#     h = b-A*a
#
#     Gp, dGp = Gdiff(x)
#
#     Fc = Ms*Gp*h
#
#     dFc1 = kron(permutedims(Gp*h),I(size(M,1)))*dMs
#     dFc2 = kron(permutedims(h),Ms)*dGp
#     dFc3 = Ms*Gp*ForwardDiff.jacobian(z->hdiff(z),x)
#
#     dFc = dFc1 + dFc2 + dFc3;
#
#     dFu = ForwardDiff.jacobian(z->assembleF(z,u,j,GravityInInertial),x)
#     return Fc, dFc, Fu, dFu
# end
# # finJ_Fc = FiniteDifferences.jacobian(m,z->FcFn(z,u,j,g)[1],x0)[1];
# # fdJ_Fc = Fcdiff(x0)[2];
# # println("\nerr_Fc = ", norm(finJ_Fc - fdJ_Fc))
# #
# # finJ_Fu = FiniteDifferences.jacobian(m,z->Fcdiff(z)[3],x0)[1]
# # fdJ_Fu = Fcdiff(x0)[4]
# # println("\nerr_Fu = ", norm(finJ_Fu - fdJ_Fu))
#
# ## fxdotdiff
# invJ = j[1].RB2.invJ;
# nJ = length(j); # also equal to number of bodies (excluding inertial frame)
# function accdiff(x::Vector{T}) where T<:Real
#     len_x = length(x)-14
#
#     M, _, _, _, _, _, M_inv, dM_inv = Mdiff(x)
#     Fc, dFc, Fu, dFu = Fcdiff(x)
#
#     dacc1 = kron(transpose(Fc+Fu),1.0I(size(M_inv,1)))*dM_inv
#     dacc2 = M_inv*(dFc + dFu)
#     dacc = dacc1 + dacc2
#
#     # dacc = zeros(T,(7*nJ,length(x)))
#     # for k=1:nJ
#     #     ct1 = 7*(k-1); ct2 = 14*k;
#     #     dacc[ct1+1:ct1+3,:] = 1/j[k].RB2.m*(dFc[ct1+1:ct1+3,:] +   dFu[ct1+1:ct1+3,:])
#     #
#     #
#     #     E = genE(x[ct2+4:ct2+7])
#     #     fdJ_E = ForwardDiff.jacobian(z->genE(z[ct2+4:ct2+7]),x)
#     #     # fdJ_trE = ForwardDiff.jacobian(z->transpose(genE(z[ct2+4:ct2+7])),x)
#     #     fdJ_trE = createCommMat(E)*fdJ_E;
#     #
#     #     dacc2_1 = kron(permutedims(invJ*E*((Fc+Fu)[ct1+4:ct1+7])),I(size(E,1)))*fdJ_trE
#     #     dacc2_2 = kron(permutedims((Fc+Fu)[ct1+4:ct1+7]),permutedims(E)*invJ)*fdJ_E
#     #     dacc2_3 = permutedims(E)*invJ*E*(dFc[ct1+4:ct1+7,:] + dFu[ct1+4:ct1+7,:])
#     #     # println()
#     #     # println("\ndacc2_1 = ", dacc2_1)
#     #     # println("\ndacc2_2 = ", dacc2_2)
#     #     # println("\ndacc2_3 = ", dacc2_3)
#     #     # sleep(1000);
#     #     dacc[ct1+4:ct1+7,:] = 1/4*(dacc2_1 + dacc2_2 + dacc2_3)
#     # end
#     return dacc
# end
# dacc = accdiff(x0);
# ##
# function fxdotdiff(x::Vector{T}) where T<:Real
#     len_x = length(x)
#     dFx = zeros(T,(len_x,len_x))
#     # dFx = Matrix{T}(undef,(len_x,len_x))
#
#     dacc = accdiff(x);
#     # println("dacc = ", dacc)
#
#     for k=1:length(j)
#         dFx[14*k+1:14*k+7,14*k+8:14*k+14] = 1.0*I(7)
#         dFx[14*k+8:14*k+14,:] = dacc[7*(k-1)+1:7*k,:]
#     end
#     return dFx
# end
#
# finJ_fx = FiniteDifferences.jacobian(m,z->fxdot(z,u,j,g),x0)[1];
# # finJ_fxM = FiniteDifferences.jacobian(m,z->fxdotM(z,u,j,g),x0)[1];
# # norm(finJ_fx - finJ_fxM)
# # norm(finJ_fx[25:28,18:21] - finJ_fxM[25:28,18:21])
#
# fdJ_fx = fxdotdiff(x0);
# println("\nerr_fx = ", norm(finJ_fx- fdJ_fx))
# ## fixing fxdot
# # function fxdotOrig(X::Vector{T},U::Matrix{S},j::Vector{Joint},GravityInInertial::Vector{Float64}) where {T<:Real,S<:Real}
# #     nJ = length(j); # Number of joints
# #
# #     unconstrF, constrF = Constraint(X, U, j, GravityInInertial)
# #
# #     # m = central_fdm(10,1)
# #     # jb1 = FiniteDifferences.jacobian(m,z -> Constraint(z,U,j,GravityInInertial)[2],X)[1]
# #     # jb2 = ForwardDiff.jacobian(z ->Constraint(z,U,j,GravityInInertial)[2],X)
# #     # println("jb_Fcerr = ", norm(jb1-jb2))
# #     # sleep(1000)
# #
# #     xdot = Vector{Union{T,S}}(undef,14*(nJ+1));
# #     xdot[1:14] = zeros(14);
# #
# #     for k=1:nJ
# #         xb = X[14*k+1:14*(k+1)]
# #
# #         unconstrainedF_rb = unconstrF[:,j[k].RB2.bodyID]
# #         Fc_rb = constrF[:,j[k].RB2.bodyID]
# #         TotalForce = unconstrainedF_rb[1:3] + Fc_rb[1:3]
# #         TotalMoment = unconstrainedF_rb[4:7] + Fc_rb[4:7]
# #
# #         # xdot = xb[8:10]
# #         E = genE(X[14*k+4:14*k+7])
# #         xdot[14*k+1:14*k+3] = X[14*k+8:14*k+10]#xb[8:10]
# #         xdot[14*k+4:14*k+7] = X[14*k+11:14*k+14]
# #         xdot[14*k+8:14*k+10] = TotalForce/j[k].RB2.m
# #         xdot[14*k+11:14*k+14] = 1/4*transpose(E)*j[k].RB2.invJ*E*TotalMoment
# #         # xdot[14*k+8:14*k+10] = TotalForce/j[k].RB2.m
# #         # xdot[14*k+11:14*k+14] = TotalMoment
# #     end
# #
# #     return xdot
# # end
# # function fxdotM(X::Vector{T},U::Matrix{S},j::Vector{Joint},GravityInInertial::Vector{Float64}) where {T<:Real,S<:Real}
# #     nJ = length(j); # Number of joints
# #
# #     unconstrF, constrF = Constraint(X, U, j, GravityInInertial)
# #
# #     # m = central_fdm(10,1)
# #     # jb1 = FiniteDifferences.jacobian(m,z -> Constraint(z,U,j,GravityInInertial)[2],X)[1]
# #     # jb2 = ForwardDiff.jacobian(z ->Constraint(z,U,j,GravityInInertial)[2],X)
# #     # println("jb_Fcerr = ", norm(jb1-jb2))
# #     # sleep(1000)
# #
# #     xdot = Vector{Union{T,S}}(undef,14*(nJ+1));
# #     xdot[1:14] = zeros(14);
# #
# #     for k=1:nJ
# #         xb = X[14*k+1:14*(k+1)]
# #
# #         unconstrainedF_rb = unconstrF[:,j[k].RB2.bodyID]
# #         Fc_rb = constrF[:,j[k].RB2.bodyID]
# #         TotalForce = unconstrainedF_rb[1:3] + Fc_rb[1:3]
# #         TotalMoment = unconstrainedF_rb[4:7] + Fc_rb[4:7]
# #         genF = [TotalForce; TotalMoment]; # generalized force acting on RB
# #         genM = genMatM(X, j[k].RB2) # generalized mass matrix of RB
# #         xdot[14*k+1:14*k+7] = X[14*k+8:14*k+14]#xb[8:10]
# #         xdot[14*k+8:14*k+14] = inv(genM)*genF;
# #         # xdot[14*k+8:14*k+10] = TotalForce/j[k].RB2.m;
# #         # Mbeta_inv = inv(genM[4:7,4:7]);
# #         # xdot[14*k+11:14*k+14] = TotalMoment;
# #     end
# #
# #     return xdot
# # end
# # xdot1 = fxdotOrig(x0,u,j,g);
# # xdot2 = fxdotM(x0,u,j,g);
# # norm(xdot1 - xdot2)
# ## M for β = [1.0;zeros(3)]
# # q = deepcopy(x0); q[18:21] = [1.0;zeros(3)]
# # function Mq(x::AbstractVector{T}) where T<:Real
# #     M = assembleM(x,j)
# #     M1 = M#[4:7,4:7]
# #     out = (M1)#\I(size(M1,1))
# #     return out
# # end
# # finJ_M = FiniteDifferences.jacobian(m,z->inv(Mq(z)), x0)[1]
# # M1 = Mq(x0);
# # fdJ_M = -(kron(permutedims(inv(M1)),inv(M1)))*ForwardDiff.jacobian(Mq, x0)
# # println("\nerr_invM = ", norm(finJ_M - fdJ_M))
#
# ## commutation matrix
# function createCommMat(M)
#     r,m = size(M);
#     K = zeros(m*r,m*r);
#     for i=1:r
#         for j=1:m
#             ei = 1.0I(r)[:,i] # ith-canonical unit vector of dimension r
#             ej = 1.0I(m)[:,j] # jth-canonical unit vector of dimension m
#             K += kron(ei*permutedims(ej), ej*permutedims(ei))
#         end
#     end
#     return K
# end
# ## checking transpose derivative
# # A = Ab_VecOfMat(x0,j)[1];
# # finJ_At = FiniteDifferences.jacobian(m,z->permutedims(Ab_VecOfMat(z,j)[1]),x0)[1]
# # fdJ_A = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x0)
# # fdJ_At = createCommMat(A)*fdJ_A;
# # println("\nerr_A = ", norm(finJ_At - fdJ_At))
#
# ## checking kroneckerProduct derivative
# # A = Ab_VecOfMat(x0,j)[1];
# # M, _, _, dM, _, _, _, _ = Mdiff(x0);
# #
# # fdJ_A = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x0)
# # finJ_AM = FiniteDifferences.jacobian(m, z->kron(Ab_VecOfMat(z,j)[1],assembleM(z,j)),x0)[1];
# function fdJ_kronAB(A,B, dA, dB)
#     n,q = size(A); p,r = size(B);
#     Iq = 1.0I(q); Ip = 1.0I(p); In = 1.0I(n); Ir = 1.0I(r);
#     Inq = Matrix{Float64}(I,n,q); Ipr = Matrix{Float64}(I,p,r)
#     Krn = createCommMat(rand(r,n));
#     Ay1 = (kron(Krn,Ip)); Ay2 = (kron(In,vec(B)));
#     # println("size(Ay1) = ", size(Ay1));
#     # println("size(Ay2) = ", size(Ay2))
#     Ay = kron(Iq, Ay1*Ay2);
#     # Ay = kron(Iq, (kron(Krn,Ip))*(kron(In,vec(B))));
#     Bx = kron((kron(Iq,Krn)*kron(vec(A),Ir)), Ip);
#     out = Ay*dA + Bx*dB
#     return out
# end
#
# function fdJ_kronAB2(A,B, dA, dB)
#     n,q = size(A); p,r = size(M);
#     Iq = 1.0I(q); Ip = 1.0I(p); In = 1.0I(n); Ir = 1.0I(r);
#     Inq = Matrix{Float64}(I,n,q); Ipr = Matrix{Float64}(I,p,r)
#     Krn = createCommMat(rand(r,n));
#     Ay1 = (kron(Krn,Ip)); Ay2 = (kron(In,vec(B)));
#     println("size(Ay1) = ", size(Ay1));
#     println("size(Ay2) = ", size(Ay2))
#     Ay = kron(Iq, Ay1*Ay2);
#     # Ay = kron(Iq, (kron(Krn,Ip))*(kron(In,vec(B))));
#     Bx = kron((kron(Iq,Krn)*kron(vec(A),Ir)), Ip);
#     out = Ay*dA + Bx*dB
#     return out
# end
# # fdJ_AM = fdJ_kronAB(A,M,fdJ_A, dM);
# # println("\nerr_AM = ", norm(finJ_AM- fdJ_AM))
#
# ## hessian of M square root
# function z1f(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x)
#     z1= kronSum(permutedims(Ms), Ms);
#     z1_inv = inv(z1);
#     dMs_t = createCommMat(Ms)*dMs;
#     dz1_1 = fdJ_kronAB(permutedims(Ms),1.0I(size(Ms,1)),dMs_t,zeros(size(dMs_t)));
#     dz1_2 = fdJ_kronAB(1.0I(size(Ms,2)),Ms,zeros(size(dMs)),dMs);
#     dz1 = dz1_1 + dz1_2;
#     dz1_inv = -(kron(permutedims(z1_inv),z1_inv))*dz1;
#     return z1_inv, dz1_inv
# end
# # finJ_z1_inv = FiniteDifferences.jacobian(m, z->z1f(z)[1],x0)[1];
# # fdJ_z1_inv = z1f(x0)[2];
# # norm(finJ_z1_inv - fdJ_z1_inv)
# # dM_f(x) = ForwardDiff.jacobian(z->assembleM(z,j),x); # function for
# # dM2 = ForwardDiff.jacobian(z->dM_f(z),x0)
# # dM2_fin = FiniteDifferences.jacobian(m, z->dM_f(z),x0)[1]
# # norm(dM2 - dM2_fin)
# function hessMs(x)
#     dM_f(y) = ForwardDiff.jacobian(z->assembleM(z,j),y);
#     dM2 = ForwardDiff.jacobian(z->dM_f(z),x); # Hessian of M
#
#     z1_inv, dz1_inv = z1f(x);
#     hessMs1 = kron(permutedims(dM), 1.0I(size(z1_inv,1)))*dz1_inv;
#     hessMs2 = kron(1.0I(size(dM,2)), z1_inv)*dM2
#     hessMs = hessMs1 + hessMs2;
#     return hessMs;
# end
# # dMs_fin(x) = FiniteDifferences.jacobian(m,z->Mdiff(z)[2],x)[1];
# # fdH_Ms = hessMs(x0);
# # finH_Ms = FiniteDifferences.jacobian(m, z->dMs_fin(z), x0)[1];
# # norm(fdH_Ms - finH_Ms)
# ## jacobian of Y1 = kron((G†h)^T, I)
# function jacY1(x)
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x0); dh = ForwardDiff.jacobian(z->hdiff(z),x0);
#     Y1 = kron(permutedims(Gp*h),1.0I(size(Ms,1)));
#     dGph = kron(permutedims(h),1.0I(size(Gp,1)))*dGp + Gp*dh;
#     dGph_t = dGph;
#     dY1 = fdJ_kronAB(permutedims((Gp*h)[:,:]),1.0I(size(Ms,1)), dGph_t, zeros(size(dMs)))
#     # dh_t = dh;
#     # dGp_t = createCommMat(Gp)*dGp;
#     # dY1 = kron(Gp, 1.0I(1))*dh_t + kron(1.0I(size(Gp,1)), permutedims(h))*dGp_t;
#     return Y1, dY1
# end
# # Gp, dGp = Gdiff(x0);
# # finJ_Y1 = FiniteDifferences.jacobian(m,z->jacY1(z)[1],x0)[1];
# # fdJ_Y1 = jacY1(x0)[2];
# # norm(finJ_Y1 - fdJ_Y1)
# ## hessian of matrix square root inverse
# function hessMsInv(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#
#     R = kron(permutedims(Ms_inv), Ms_inv);
#     dMs_inv_t = createCommMat(Ms_inv)*dMs_inv;
#     dR = fdJ_kronAB(permutedims(Ms_inv), Ms_inv, dMs_inv_t, dMs_inv);
#
#     dMs2 = hessMs(x);
#
#     out1 = kron(permutedims(dMs),1.0I(size(R,1)))*dR;
#     out2 = kron(1.0I(size(dMs,2)),R)*dMs2;
#     hessMsInv = -(out1 + out2);
#
#     return hessMsInv;
# end
# # dMsInv_fin(x) = FiniteDifferences.jacobian(m, z->Mdiff(z)[3],x)[1];
# # finH_MsInv = FiniteDifferences.jacobian(m, z->dMsInv_fin(z), x0)[1];
# # fdH_MsInv = hessMsInv(x0);
# # norm(fdH_MsInv - finH_MsInv)
# ## hessian of A
# function jacA(x)
#     dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
#     return dA
# end
# function hessA(x)
#     dA(y) = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],y);
#     cfg1 = ForwardDiff.JacobianConfig(dA, x, ForwardDiff.Chunk{1}());
#     dA2 = ForwardDiff.jacobian(dA,x, cfg1);
#     return dA2
# end
# # dA_fin(y) = FiniteDifferences.jacobian(m, z->Ab_VecOfMat(z,j)[1], y)[1];
# # finH_A = FiniteDifferences.jacobian(m, z->dA_fin(z),x0)[1];
# # fdH_A = hessA(x0);
# # norm(fdH_A- finH_A)
#
# ## hessian of G = (AM^{-1/2})
# function hessG(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     A = Ab_VecOfMat(x,j)[1]; dA = jacA(x);
#     Q1 = kron(permutedims(Ms_inv), 1.0I(size(A,1)));
#     dMs_inv_t = createCommMat(Ms_inv)*dMs_inv;
#     dQ1 = fdJ_kronAB(permutedims(Ms_inv), 1.0I(size(A,1)),dMs_inv_t,zeros(size(A,1)^2,length(x0)));
#
#     fdH_A = hessA(x);
#
#     fdH_MsInv = hessMsInv(x);
#
#     Q2 = kron(1.0I(size(Ms_inv,2)), A);
#     dQ2 = fdJ_kronAB(1.0I(size(Ms_inv,2)), A, zeros(size(dMs_inv)), dA);
#
#
#     out1 = kron(permutedims(dA),1.0I(size(Q1,1)))*dQ1;
#     out2 = kron(1.0I(size(dA,2)), Q1)*fdH_A;
#     out3 = kron(permutedims(dMs_inv), 1.0I(size(Q2,1)))*dQ2;
#     out4 = kron(1.0I(size(dMs_inv,2)), Q2)*fdH_MsInv;
#
#     hessG = out1 + out2 + out3 + out4;
#     return hessG
# end
# function G(x)
#     A = Ab_VecOfMat(x,j)[1];
#     M = assembleM(x,j);
#     G = A*M^(-1/2);
#     return G
# end
# function jacG(x)
#     A = Ab_VecOfMat(x,j)[1];
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     G = A*M^(-1/2);
#     dA = jacA(x);
#     dG1 = kron(permutedims(Ms_inv),1.0I(size(A,1)))*dA;
#     dG2 = kron(1.0I(size(Ms_inv,2)),A)*dMs_inv;
#     dG = dG1 + dG2;
#     return dG;
# end
# # dG_fin(y) = FiniteDifferences.jacobian(m, z->G(z),y)[1];
# # finH_G = FiniteDifferences.jacobian(m, z->dG_fin(z), x0)[1];
# # fdH_G = hessG(x0);
# # norm(fdH_G - finH_G)
# ## hess of H = G*(G^T)
# function hessH(x)
#     A = Ab_VecOfMat(x,j)[1];
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     G = A*M^(-1/2);
#
#     dG = jacG(x);
#     dG_t = createCommMat(G)*dG;
#
#     dH = kron(G,1.0I(size(G,1)))*dG + kron(1.0I(size(G,1)),G)*dG_t;
#
#     P1 = kron(G,1.0I(size(G,1)));
#     dP1 = fdJ_kronAB(G, 1.0I(size(G,1)), dG, zeros(size(G,1)^2,length(x)));
#
#     P2 = kron(1.0I(size(G,1)), G);
#     dP2 = fdJ_kronAB(1.0I(size(G,1)), G, zeros(size(G,1)^2,length(x)), dG);
#
#     dG2 = hessG(x);
#     dG2_t = kron(1.0I(size(dG,2)), createCommMat(G))*dG2;
#
#     out1 = kron(permutedims(dG), 1.0I(size(P1,1)))*dP1;
#     out2 = kron(1.0I(size(dG,2)), P1)*dG2;
#     out3 = kron(permutedims(dG_t), 1.0I(size(P2,1)))*dP2;
#     out4 = kron(1.0I(size(dG_t,2)), P2)*dG2_t;
#
#     fdH_H = out1 + out2 + out3 + out4;
#     return dH, fdH_H
# end
# function fnH(x)
#     A = Ab_VecOfMat(x,j)[1];
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     G = A*M^(-1/2);
#     H = G*permutedims(G);
#     return H
# end
# # fdJ_H, fdH_H = hessH(x0);
# # dH_fin(y) = FiniteDifferences.jacobian(m, z->fnH(z), y)[1];
# # finH_H = FiniteDifferences.jacobian(m, z->dH_fin(z), x0)[1];
# # norm(dH_fin(x0) - fdJ_H)
# # norm(finH_H - fdH_H)
# ## hess of G†
# function hessGp(x)
#     H = fnH(x); H_inv = inv(H);
#     dH, dH2 = hessH(x);
#     Gnum = G(x);
#     dG = jacG(x); dG2 = hessG(x);
#     dG_t = createCommMat(Gnum)*dG;
#     dG2_t = kron(1.0I(size(dG,2)), createCommMat(Gnum))*dG2;
#
#     dH_inv = -kron(permutedims(H_inv), H_inv)*dH;
#     dH_inv_t = createCommMat(H_inv)*dH_inv;
#
#     W = kron(permutedims(H_inv), H_inv);
#     dW = fdJ_kronAB(permutedims(H_inv), H_inv, dH_inv_t, dH_inv);
#     h1 = kron(permutedims(dH),1.0I(size(W,1)))*dW;
#     h2 = kron(1.0I(size(dH,2)),W)*dH2;
#     dH_inv2 = -(h1 + h2);
#
#
#     V1 = kron(permutedims(H_inv), 1.0I(size(Gnum,2)));
#     dV1 = fdJ_kronAB(permutedims(H_inv), 1.0I(size(Gnum,2)), dH_inv_t, zeros(size(Gnum,2)^2, length(x0)));
#
#     V2 = kron(1.0I(size(H_inv,2)), permutedims(Gnum));
#     dV2 = fdJ_kronAB(1.0I(size(H_inv,2)),permutedims(Gnum), zeros(size(dH_inv_t)), dG_t);
#
#     dGp = V1*dG_t + V2*dH_inv;
#
#     out1 = kron(permutedims(dG_t), 1.0I(size(V1,1)))*dV1;
#     out2 = kron(1.0I(size(dG_t,2)), V1)*dG2_t;
#     out3 = kron(permutedims(dH_inv), 1.0I(size(V2,1)))*dV2;
#     out4 = kron(1.0I(size(dH_inv,2)), V2)*dH_inv2;
#
#     dGp2 = out1 + out2 + out3 + out4;
#
#     return dGp2;
# end
# # fdH_Gp = hessGp(x0);
# # dGp_fin(y) = FiniteDifferences.jacobian(m,z->pinv(G(z)),y)[1];
# # finH_Gp = FiniteDifferences.jacobian(m, z->dGp_fin(z), x0)[1];
# # norm(finH_Gp - fdH_Gp)
#
# ## hessian of h
# # dh(y) = ForwardDiff.jacobian(z->hdiff(z),y);
# # dh_fin(y) = FiniteDifferences.jacobian(m, z->hdiff(z), y)[1];
# # norm(dh(x0) - dh_fin(x0))
# # cfg1 = ForwardDiff.JacobianConfig(dh, x0, ForwardDiff.Chunk{1}());
# # fdH_h = ForwardDiff.jacobian(dh,x0,cfg1);
# # finH_h = FiniteDifferences.jacobian(m,z->dh_fin(z),x0)[1];
# # norm(fdH_h - finH_h)
#
# ## jacobian of T1 = ((G†h)^T ⊗ I)*dMs
# function jacT1(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x);
#     jac_h = ForwardDiff.jacobian(z->hdiff(z),x);
#
#     D1 = kron(permutedims(Gp*h), 1.0I(size(Ms,1)));
#     dGph = kron(permutedims(h), 1.0I(size(Gp,1)))*dGp + Gp*jac_h;
#     dGph_t = dGph;
#     dD1 = fdJ_kronAB(permutedims(Gp*h), 1.0I(size(Ms,1)), dGph_t, zeros(size(dMs)));
#
#     dMs2 = hessMs(x);
#     dMs_t = createCommMat(Ms)*dMs;
#
#     dT1_1 = kron(permutedims(dMs), 1.0I(size(D1,1)))*dD1;
#     dT1_2 = kron(1.0I(size(dMs,2)), D1)*dMs2;
#
#     dT1 = dT1_1 + dT1_2;
#     return dT1
# end
# function fnT1(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x);
#     D1 = kron(permutedims(Gp*h), 1.0I(size(Ms,1)));
#
#     T1 = D1*dMs;
#     return T1
# end
# # fdJ_T1 = jacT1(x0);
# # finJ_T1 = FiniteDifferences.jacobian(m, z->fnT1(z), x0)[1];
# # norm(fdJ_T1 - finJ_T1)
#
# ## jacobian of T2
# function fnT2(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x);
#     jac_h = ForwardDiff.jacobian(z->hdiff(z),x);
#
#     D2 = kron(permutedims(h), Ms);
#     T2 = D2*dGp;
#     return T2
# end
#
# function jacT2(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x);
#     jac_h = ForwardDiff.jacobian(z->hdiff(z),x);
#
#     D2 = kron(permutedims(h), Ms);
#     dD2 = fdJ_kronAB(permutedims(h), Ms, jac_h, dMs);
#
#     dGp2 = hessGp(x);
#
#     dT2_1 = kron(permutedims(dGp),1.0I(size(D2,1)))*dD2;
#     dT2_2 = kron(1.0I(size(dGp,2)), D2)*dGp2;
#     dT2 = dT2_1 + dT2_2;
#     return dT2
# end
# # fdJ_T2 = jacT2(x0);
# # finJ_T2 = FiniteDifferences.jacobian(m, z->fnT2(z), x0)[1];
# # norm(fdJ_T2 - finJ_T2)
# ## jacobian of T3
# function fnT3(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x);
#     jac_h = ForwardDiff.jacobian(z->hdiff(z),x);
#
#     T3 = Ms*Gp*jac_h;
#     return T3
# end
# function jacT3(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h = hdiff(x);
#     # jac_h = ForwardDiff.jacobian(z->hdiff(z),x);
#     dh(y) = ForwardDiff.jacobian(z->hdiff(z),y);
#     jac_h = dh(x);
#     cfg1 = ForwardDiff.JacobianConfig(dh, x0, ForwardDiff.Chunk{1}());
#     dh2 = ForwardDiff.jacobian(dh,x,cfg1);
#
#     dT3_1 = kron(permutedims(Gp*jac_h), 1.0I(size(Ms,1)))*dMs;
#     dT3_2 = kron(permutedims(jac_h), Ms)*dGp;
#     dT3_3 = kron(1.0I(size(jac_h,2)), Ms*Gp)*dh2;
#     dT3 = dT3_1 + dT3_2 + dT3_3;
#
#     return dT3
# end
# function hessFc(x)
#     # second partial derivative of fc wrt x, i.e., d/dx(dfc/dx)
#     out = jacT1(x) + jacT2(x) + jacT3(x)
#     return out
# end
#
# # fdJ_T3 = jacT3(x0);
# # finJ_T3 = FiniteDifferences.jacobian(m, z->fnT3(z), x0)[1];
# # norm(fdJ_T3 - finJ_T3)
# # fdH_Fc = fdJ_T1 + fdJ_T2 + fdJ_T3;#hessFc(x0);
# # finH_Fc = FiniteDifferences.jacobian(m, z->Fcdiff(z)[2], x0)[1];
# # norm(fdH_Fc - finH_Fc)
#
# ## assembleFdiff
# function assembleFdiff(x, u, j::Vector{Joint}, GravityInInertial:: Vector{Float64})
#     maxBodyID = j[end].RB2.bodyID
#     F = Vector(undef,7*(maxBodyID-1))
#
#     for i=1:length(j)
#         k = j[i].RB2.bodyID - 1
#         F[7*(k-1)+1:7*k] =
#         genExtF(x,j[i].RB2,u[:,k+1],GravityInInertial)
#
#         if j[i].type == "Spring"
#             F_b1, F_b2 = ForceConSpr(x,j[i]); # Forces acting on the bodies connected by the spring. Note that a spring does not constrain 2 bodies any more than they previously were, i.e., degrees of freedom remain the same.
#
#             b1_id = j[i].RB1.bodyID -1; b2_id = j[i].RB2.bodyID - 1;
#             if j[i].RB1.m == 0.0 # First body is the inertial frame
#                 F[7*(b2_id-1)+1:7*b2_id] += F_b2;
#             else # First body is not the inertial frame
#                 F[7*(b1_id-1)+1:7*b1_id] += F_b1;
#                 F[7*(b2_id-1)+1:7*b2_id] += F_b2;
#             end
#         end
#     end
#     return F
# end
# ## hessian of Fu
# # dFu(y) = ForwardDiff.jacobian(z->assembleF(z,u,j,g), y);
# # dFu_fin(y) = FiniteDifferences.jacobian(m, z->assembleF(z,u,j,g), y)[1];
# # # norm(dFu(x0) - dFu_fin(x0))
# # fdH_Fu = ForwardDiff.jacobian(z->dFu(z),x0);
# # finH_Fu = FiniteDifferences.jacobian(m, z->dFu_fin(z),x0)[1];
# # norm(fdH_Fu - finH_Fu)
# ##
# ##
# ## d/du(dh/dx)
# ## d/du(dh/dx)
# function hdiff_u(x,u)
#     Fu = assembleFdiff(x,u,j,g);
#     A, b = Ab_VecOfMat(x,j);
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     h = b - A*inv(M)*Fu;
#
#     dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
#     db = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[2],x);
#     dFu = dFu_x_fn(x,u);
#     dFu_u = ForwardDiff.jacobian(z->assembleFdiff(x,z,j,g),u);
#
#     # dh_x
#     term1 = kron(permutedims(M_inv*Fu), 1.0I(size(A,1)))*dA;
#     term2 = kron(permutedims(Fu), A)*dM_inv;
#     term3 = A*M_inv*dFu;
#     dh_x = db - (term1 + term2 + term3);
#
#     # dh_u
#     dh_u = -A*M_inv*dFu_u;
#
#     # dh_u_x
#     dFu_u_x = ForwardDiff.jacobian(y->dFu_x_fn(x,y),u);
#     E1 = kron(permutedims(M_inv*Fu), 1.0I(size(A,1)));
#     dE1_u = fdJ_kronAB(permutedims(M_inv*Fu), 1.0I(size(A,1)), M_inv*dFu_u, zeros(size(A,1)^2, length(u)));
#
#
#     E2 = kron(permutedims(Fu), A);
#     dE2_u = fdJ_kronAB(permutedims(Fu), A, dFu_u, zeros(length(A), length(u)));
#
#
#     term1 = kron(permutedims(dA),1.0I(size(E1,1)))*dE1_u;
#     term2 = kron(permutedims(dM_inv), 1.0I(size(E2,1)))*dE2_u;
#     term3 = kron(1.0I(size(dFu,2)), A*M_inv)*dFu_u_x;
#
#     dh_u_x = -(term1 + term2 + term3);
#
#     # dh_x_u
#     dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
#     dFu_x_u = ForwardDiff.jacobian(y->dFu_u_fn(y,u),x);
#
#     out1 = kron(permutedims(M_inv*dFu_u), 1.0I(size(A,1)))*dA;
#     out2 = kron(permutedims(dFu_u), A)*dM_inv;
#     out3 = kron(1.0I(length(u)), A*M_inv)*dFu_x_u;
#     dh_x_u = -(out1 + out2 + out3);
#     return h, dh_u, dh_u_x, dh_x_u
# end
# fdH_h_xu = hdiff_u(x0,u)[4];
# # fdH_h_xu = hdiff_u(x0,u)[2];
# # fin_h_x(y, uVar) = FiniteDifferences.jacobian(m, z->hdiff_u(z,uVar)[1], y)[1];
# # fin_h_xu = FiniteDifferences.jacobian(m, z->fin_h_x(x0,z),u)[1];
# # dh_u_temp = ForwardDiff.jacobian(z->hdiff_u(x0,z)[2],u);
# # norm(fdH_h_xu - fin_h_xu)
# ## d/du(dFu/dx)
# function dFu_x_fn(x,u)
#     # x = z[1:28]; u = reshape(z[29:end],(6,2))
#     out = ForwardDiff.jacobian(y->assembleFdiff(y,u,j,g),x)
#     return out
# end
# function dFu_u_fn(x,u)
#     out = ForwardDiff.jacobian(y->assembleFdiff(x,y,j,g),u);
#     return out
# end
# # dFu_x(x0,u)
# # dFu_x_u = ForwardDiff.jacobian(y->dFu_x(x0,y),u);
# # fin_dFu_x(y,uVar) = FiniteDifferences.jacobian(m, z->assembleF(z,uVar,j,g),y)[1];
# # fin_dFu_x_u = FiniteDifferences.jacobian(m, z->fin_dFu_x(x0,z), u)[1];
# # norm(dFu_x(x0,u) - fin_dFu_x(x0,u))
# # norm(dFu_x_u - fin_dFu_x_u)
# ## d/du(dFu/du)
# function dFu_uu(x,u)
#     dFu_u(y) = ForwardDiff.jacobian(z->assembleFdiff(x,z,j,g),y);
#     out = ForwardDiff.jacobian(z->dFu_u(z), u);
#     return out
# end
# # hessFu_uu = dFu_uu(x0,u);
# # fin_dFu_u(y,uVar) = FiniteDifferences.jacobian(m, z->assembleF(y,z,j,g),uVar)[1];
# # fin_dFu_uu = FiniteDifferences.jacobian(m, z->fin_dFu_u(x0,z), u)[1];
# # norm(hessFu_uu-fin_dFu_uu)
# ## d/dx(dFu/dx)
# function dFu_xx(x,u)
#     dFu_x(y) = ForwardDiff.jacobian(z->assembleFdiff(z,u,j,g),y);
#     out = ForwardDiff.jacobian(z->dFu_x(z), x);
#     return out
# end
# # hessFu_xx = dFu_xx(x0,u);
# # fin_dFu_x(y) = FiniteDifferences.jacobian(m, z->assembleF(z,u,j,g), y)[1];
# # fin_dFu_xx = FiniteDifferences.jacobian(m, z->fin_dFu_x(z), x0)[1];
# # norm(hessFu_xx - fin_dFu_xx)
# ## d/du(dfc/dx)
# function cross_fc_ux(x,u)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h, dh_u, dh_ux, dh_xu = hdiff_u(x,u);
#     dh(y) = ForwardDiff.jacobian(z->hdiff(z),y);
#     jac_h = dh(x);
#
#     # first term
#     D1 = kron(permutedims(Gp*h), 1.0I(size(Ms,1)));
#     dGph_u = Gp*dh_u;
#     dGph_u_t = dGph_u;
#     dD1_u = fdJ_kronAB(permutedims(Gp*h), 1.0I(size(Ms,1)), dGph_u_t, zeros(size(dMs,1),length(u)));
#     # dMs_t = createCommMat(Ms)*dMs;
#     out1 = kron(permutedims(dMs), 1.0I(size(D1,1)))*dD1_u;
#
#     # second term
#     D2 = kron(permutedims(h), Ms)*dGp;
#     dD2_u = fdJ_kronAB(permutedims(h), Ms, dh_u, zeros(size(dMs,1), length(u)));
#     out2 = kron(permutedims(dGp), 1.0I(size(D2,1)))*dD2_u;
#
#     # third term
#     out3 = kron(1.0I(size(jac_h,2)), Ms*Gp)*dh_ux;
#
#     d2Fc_xu = out1 + out2 + out3;
#     return d2Fc_xu;
# end
#
# function cross_fc_xu(x,u)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#     h, dh_u, dh_ux, dh_xu = hdiff_u(x,u);
#
#     out1 = kron(permutedims(Gp*dh_u), 1.0I(size(Ms,1)))*dMs;
#     out2 = kron(permutedims(dh_u), Ms)*dGp;
#     out3 = kron(1.0I(length(u)), Ms*Gp)*dh_xu;
#     out = out1 + out2 + out3;
#     return out
# end
# fdH_Fc_xu = cross_fc_xu(x0,u);
# # fin_dFc_x(y,uVar) = FiniteDifferences.jacobian(m, z->fnFc(z,uVar), y)[1];
# fin_dFc_u(xVar,uVar) = FiniteDifferences.jacobian(m, z->fnFc(xVar,z), uVar)[1];
# fin_dFc_xu = FiniteDifferences.jacobian(m, z->fin_dFc_u(z,u),x0)[1];
# # fin_dFc_ux = FiniteDifferences.jacobian(m, z->fin_dFc_x(x0,z), u)[1];
# norm(fdH_Fc_xu - fin_dFc_xu)
# ## d/du(dFc/du)
# function hessFc_uu(x,u)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     A,b = Ab_VecOfMat(x,j);
#     Gp, dGp = Gdiff(x);
#     h, dh_u, dh_ux, dh_xu = hdiff_u(x,u);
#     hessFu_uu = dFu_uu(x,u);
#
#     dh_uu = -kron(1.0I(length(u)), A*M_inv)*hessFu_uu;
#
#     out = kron(1.0I(size(dh_u,2)), Ms*Gp)*dh_uu;
#     return out
# end
# # dFc_uu = hessFc_uu(x0,u);
# function fnFc(x,u)
#     Fu = assembleF(x,u,j,GravityInInertial)
#     A, b = Ab_VecOfMat(x,j)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#
#     a = M_inv*Fu
#     h = b-A*a
#
#     Gp, dGp = Gdiff(x)
#
#     Fc = Ms*Gp*h;
#     return Fc
# end
# # fin_dFc_u(y) = FiniteDifferences.jacobian(m, z->fnFc(x0,z), y)[1];
# # fin_dFc_uu = FiniteDifferences.jacobian(m, z->fin_dFc_u(z), u)[1];
# # norm(fin_dFc_uu - dFc_uu)
# ## beta component of acceleration only
# function hessqdd(x)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Fc, dFc, Fu, dFu = Fcdiff(x)
#
#     # Hessian of M and M_inv
#     dM_f(y) = ForwardDiff.jacobian(z->assembleM(z,j),y);
#     dM2 = ForwardDiff.jacobian(z->dM_f(z),x); # Hessian of M
#
#     Q = kron(permutedims(M_inv), M_inv);
#     dM_inv_t = createCommMat(M_inv)*dM_inv;
#     dQ = fdJ_kronAB(permutedims(M_inv), M_inv, dM_inv_t, dM_inv);
#     out1 = kron(permutedims(dM), 1.0I(size(Q,1)))*dQ;
#     out2 = kron(1.0I(size(dM,2)), Q)*dM2;
#     dM_inv2 = -(out1 + out2);
#
#     T1 = kron(permutedims(Fc + Fu), 1.0I(size(M_inv,1)));
#     dT1 = fdJ_kronAB(permutedims(Fc + Fu), 1.0I(size(M_inv,1)), (dFc + dFu), zeros(size(dM_inv)));
#     T2 = dFc + dFu;
#     dT2 = hessFc(x) + dFu_xx(x,u);
#
#     out1 = kron(permutedims(dM_inv), 1.0I(size(T1,1)))*dT1;
#     out2 = kron(1.0I(size(dM_inv,2)), T1)*dM_inv2;
#     out3 = kron(permutedims(T2), 1.0I(size(M_inv,1)))*dM_inv;
#     out4 = kron(1.0I(size(T2,2)), M_inv)*dT2;
#
#     out = out1 + out2 + out3 + out4;
#     return out
# end
# # fin_dM_inv(y) = FiniteDifferences.jacobian(m, z->inv(assembleM(z,j)), y)[1];
# # fin_dM_inv2 = FiniteDifferences.jacobian(m, z->fin_dM_inv(z), x0)[1];
# # hess_M_inv = hessqdd(x0);
# # norm(fin_dM_inv2 - hess_M_inv)
# function cross_qdd_ux(x,u)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Fc, dFc, Fu, dFu = Fcdiff(x)
#     Gp, dGp = Gdiff(x);
#
#     T1 = kron(permutedims(Fc + Fu), 1.0I(size(M_inv,1)));
#     dFu_u = ForwardDiff.jacobian(z->assembleFdiff(x,z,j,g),u);
#     h, dh_u, dh_ux, dh_xu = hdiff_u(x,u);
#     dFc_u = Ms*Gp*dh_u;
#     dT1_u = fdJ_kronAB(permutedims(Fc + Fu), 1.0I(size(M_inv,1)), dFc_u + dFu_u, zeros(size(M_inv,1)^2, length(u)));
#
#     T2 = dFc + dFu;
#     dFu_x_u = ForwardDiff.jacobian(y->dFu_x_fn(x,y),u);
#     dT2_u = cross_fc(x,u) + dFu_x_u;
#
#     out1 = kron(permutedims(dM_inv), 1.0I(size(T1,1)))*dT1_u;
#     out2 = kron(1.0I(length(x)), M_inv)*dT2_u;
#     out = out1+out2;
#     return out
# end
# function cross_qdd_xu(x,u)
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     h, dh_u, dh_ux, dh_xu = hdiff_u(x,u);
#     Gp, dGp = Gdiff(x);
#
#     dFc_u = Ms*Gp*dh_u;
#     dFu_u = ForwardDiff.jacobian(z->assembleFdiff(x,z,j,g),u);
#
#     dFu_xu = ForwardDiff.jacobian(y->dFu_u_fn(y,u),x);
#     dFc_xu = cross_fc_xu(x,u);
#
#     out1 = kron(permutedims(dFc_u + dFu_u), 1.0I(size(M_inv,1)))*dM_inv;
#     out2 = kron(1.0I(length(u)), M_inv)*(dFc_xu + dFu_xu);
#     out = out1 + out2;
#
#     return out
# end
# function hessqdd_u(x,u)
#     M = assembleM(x,j);
#     M_inv = inv(M);
#     d2Fc_uu = hessFc_uu(x,u);
#     d2Fu_uu = dFu_uu(x,u);
#     out = kron(1.0I(length(u)), M_inv)*(d2Fc_uu + d2Fu_uu);
#     return out
# end
# function qFn(x,u)
#     xdot = fxdot(x,u,j,g); nB = length(j)+1;
#     y = Vector(undef,7*nB); # qdd
#     for i=1:nB
#       y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
#     end
#     deleteat!(y,1:7); # inertial frame has no acceleration
#     return y
# end
# function dqFn(x,u)
#     q = qFn(x,u);
#     M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff(x);
#     Gp, dGp = Gdiff(x);
#
#     dq_x = accdiff(x);
#
#     h, dh_u, dh_ux, dh_xu = hdiff_u(x,u);
#     dFc_u = Ms*Gp*dh_u;
#     dFu_u = ForwardDiff.jacobian(z->assembleFdiff(x,z,j,g),u);
#
#     dq_u = M_inv*(dFc_u + dFu_u);
#
#     return dq_x, dq_u
# end
# # fdJ_qx, fdJ_qu = dqFn(x0,u);
# # fdJ_q = [fdJ_qx fdJ_qu];
#
# finJ_q(xVar, uVar) = FiniteDifferences.jacobian(m, z->qFn(z,uVar), xVar)[1];
# finJ_q_u(xVar, uVar) = FiniteDifferences.jacobian(m, z->qFn(xVar,z), uVar)[1];
#
# # finH_q = FiniteDifferences.jacobian(m, z->finJ_q(z,u), x0)[1];
# # fdH_q = hessqdd(x0);
# # norm(fdH_q - finH_q)
#
# # finH_q_ux = FiniteDifferences.jacobian(m, z->finJ_q(x0,z), u)[1];
# # fdH_q_ux = cross_qdd_ux(x0,u);
# # norm(finH_q_ux - fdH_q_ux)
#
# # finH_q_xu = FiniteDifferences.jacobian(m, z->finJ_q_u(z,u), x0)[1];
# # fdH_q_xu = cross_qdd_xu(x0,u);
# # norm(finH_q_xu - fdH_q_xu)
#
# # finH_q_uu = FiniteDifferences.jacobian(m, z->finJ_q_u(x0,z), u)[1];
# # fdH_q_uu = hessqdd_u(x0,u);
# # norm(finH_q_uu - fdH_q_uu)
# ## other constraints - quatNorm and JointLocation
# function quatNormConstr(z)
#     nB = length(j)+1;
#     x = z[1:14*nB]; u = reshape(z[14*nB+1:end],(6,nB));
#     q_constr = Vector(undef, nB)
#     for i=1:nB
#         q_constr[i] = norm(x[14*(i-1)+4:14*(i-1)+7])^2 - 1;
#     end
#     return q_constr
# end
# function quatNormConstr_diff(z)
#     f(y) = sum(quatNormConstr(y).*quatNormConstr(y));
#     df(x) = ForwardDiff.gradient(y->f(y), x);
#     grad_f = df(z);
#     d2f = ForwardDiff.jacobian(y->df(y), z);
#     return grad_f, d2f
# end
# # finG_qNorm(z) = FiniteDifferences.grad(m, y->sum(quatNormConstr(y).*quatNormConstr(y)), z)[1];
# # finH_qNorm = FiniteDifferences.jacobian(m, finG_qNorm, optVar)[1];
# # fdG_qNorm, fdH_qNorm = quatNormConstr_diff(optVar);
# # norm(fdG_qNorm - finG_qNorm(optVar))
# # norm(fdH_qNorm - finH_qNorm)
# ## jointLocConstr
# function jointLocConstr(z)
#     nB = length(j)+1;
#     x = z[1:14*nB]; u = reshape(z[14*nB+1:end],(6,nB));
#     jLoc_constr = Matrix(undef,3,length(j))
#     for i=1:length(j)
#         b1_id = j[i].RB1.bodyID; b2_id = j[i].RB2.bodyID;
#         x1 = x[(b1_id-1)*14+1:b1_id*14]
#         x2 = x[(b2_id-1)*14+1:b2_id*14]
#         j_b1 = x1[1:3] + transpose(quat2dcm(x1[4:7]))*j[i].pos1;
#         j_b2 = x2[1:3] + transpose(quat2dcm(x2[4:7]))*j[i].pos2;
#         jLoc_constr[:,i] = j_b1 - j_b2;
#     end
#     return jLoc_constr
# end
# function jointLocConstr_diff(z)
#     f(y) = sum(jointLocConstr(y).*jointLocConstr(y));
#     df(x) = ForwardDiff.gradient(y->f(y), x);
#     grad_f = df(z);
#     d2f = ForwardDiff.jacobian(y->df(y), z);
#     return grad_f, d2f
# end
# # finG_jointLoc(z) = FiniteDifferences.grad(m, y->sum(jointLocConstr(y).*jointLocConstr(y)), z)[1];
# # finH_jointLoc = FiniteDifferences.jacobian(m, finG_jointLoc, optVar)[1];
# # fdG_jointLoc, fdH_jointLoc = jointLocConstr_diff(optVar);
# # norm(fdG_jointLoc - finG_jointLoc(optVar))
# # norm(fdH_jointLoc - finH_jointLoc)
# ## costFn
# function dynConstr(z)
#     nB = length(j)+1;
#     x = z[1:14*nB]; u = reshape(z[14*nB+1:end],(6,nB));
#     y = qFn(x,u);
#     out = sum(y.*y);
#     return out;
# end
# function dynConstr_jac(z)
#     jac_qdd = ForwardDiff.gradient(x->dynConstr(x), z);
#     return jac_qdd
# end
# function costFn(z)
#     nB = length(j)+1;
#     x = z[1:14*nB]; u = reshape(z[14*nB+1:end-1],(6,nB));
#
#     # Dynamics constraint
#     y = qFn(x,u);
#     # Quaternion norm constraint
#     q_constr = quatNormConstr(z[1:end-1]);
#     # JointLocation constraint  (translation constraint)
#     jLoc_constr = jointLocConstr(z[1:end-1]);
#
#     constr = [y;q_constr;jLoc_constr;]
#     J = sum(constr.*constr) - z[end];
#
#     return J
# end
# function gradCost(z)
#     nB = length(j)+1;
#     x = z[1:14*nB]; u = reshape(z[14*nB+1:end-1],(6,nB));
#     xu_var = z[1:end-1]
#     # Dynamics constraint
#     qdd = qFn(x,u);
#     # Quaternion norm constraint
#     q_constr = quatNormConstr(xu_var);
#     # JointLocation constraint  (translation constraint)
#     jLoc_constr = jointLocConstr(xu_var);
#
#     dqdd_constr = dynConstr_jac(xu_var); # Dyn
#     dq_constr,_ = quatNormConstr_diff(xu_var); # Quat
#     dj_constr, _ = jointLocConstr_diff(xu_var); # JointLoc
#
#     jac_constr = dqdd_constr + dq_constr + dj_constr;
#     out = [jac_constr;-1.0];
#     return out
# end
#
# function hessCost(z)
#     zLen = length(z);
#     nB = length(j)+1;
#     x = z[1:14*nB]; u = reshape(z[14*nB+1:end-1],(6,nB));
#     y = qFn(x,u);
#
#     dy_x, dy_u = dqFn(x,u);
#     dy = [dy_x dy_u];
#
#     d2y_xx = hessqdd(x);
#     d2y_xu = cross_qdd_xu(x,u);
#     d2y_ux = cross_qdd_ux(x,u);
#     d2y_uu = hessqdd_u(x,u);
#
#     d2y = [d2y_xx d2y_ux; d2y_xu d2y_uu];
#
#     out1 = kron(1.0I(zLen-1), permutedims(y))*d2y;
#     out2 = permutedims(dy)*dy;
#     hess_dyn = 2*(out1 + out2); # hessian of dynamics cost = sum(y.*y)
#
#     _, hess_qNorm = quatNormConstr_diff(z[1:end-1]); # hessian of quatNormConstr cost
#     _, hess_jLoc = jointLocConstr_diff(z[1:end-1]); # hessian of jointLocConstr cost
#
#     hess_constr = hess_dyn +  hess_qNorm + hess_jLoc;
#
#     out = zeros(zLen, zLen);
#     out[1:zLen-1, 1:zLen - 1] = hess_constr;
#
#     return out
# end
# optVar = [x0;u[:];0.0];
# # J = costFn(optVar);
# # Gradient comparison
# finG_J(y) = FiniteDifferences.grad(m, z->costFn(z), y)[1];
# fdG_J = gradCost(optVar);
# norm(finG_J(optVar) - fdG_J)
# # Cost Hessian comparison
#
# finH_J = FiniteDifferences.jacobian(m, finG_J, optVar)[1];
# eval_h = hessCost(optVar);
# norm(finH_J - eval_h)

## kronecker product
# using Kronecker
# using LazyArrays
# using BenchmarkTools
# # using Einsum
# # using TensorOperations
# sz = 20;
# tempA = rand(sz,sz);
# tempB = (rand(sz,sz));
# dAtemp = rand(sz^2, 10);
# dBtemp = rand(sz^2, 10);
# function kronProd(A,B)
#     A = A[:,:]; B = B[:,:];
#     # @einsum K[i, j] = A[i÷p,j÷q]*B[(i-1)%p + 1 ,(j-1)%q + 1];
#     K = ApplyArray(kron,A[:,:],B[:,:]);
#     # out = Matrix{Float64}(undef,size(K)...); copyto!(out,K);
#     # return out
#     return K
# end
# function kronProdFor(A,B)
#     A = A[:,:]; B = B[:,:];
#     rA, cA = size(A); rB, cB = size(B);
#     # K = A[1,1]*B;
#     K = Matrix(undef, rA*rB, cA*cB);
#     for i=1:rA
#         for j=1:cA
#             K[(i-1)*rB + 1: (i)*rB, (j-1)*cB + 1: j*cB] = A[i,j]*B;
#         end
#     end
#     return K
# end
# k1 = kronProd(tempA,tempB);
# k2 = kron(tempA,tempB);
# println()
# println("err = ", norm(k1-k2))
# ##
# @btime kronProd(tempA,tempB);
# @btime kron(tempA,tempB);
# ##
#
# function kronSum(A,B)
#   # kronecker sum of 2 arrays
#   out = kronProd(A, Matrix(I,(size(B)))) + kronProd(Matrix(I,(size(A))),B)
#   return out
# end
#
# function createCommMat(M)
#     # commutation matrix
#     r,m = size(M);
#     K = zeros(m*r,m*r);
#     for i=1:r
#         for j=1:m
#             ei = 1.0I(r)[:,i] # ith-canonical unit vector of dimension r
#             ej = 1.0I(m)[:,j] # jth-canonical unit vector of dimension m
#             K += kronProd(ei*permutedims(ej), ej*permutedims(ei))
#         end
#     end
#     return K
# end
#
# function fdJ_kronAB2(A,B, dA, dB)
#     n,q = size(A); p,r = size(B);
#     Iq = 1.0I(q); Ip = 1.0I(p); In = 1.0I(n); Ir = 1.0I(r);
#     Inq = Matrix{Float64}(I,n,q); Ipr = Matrix{Float64}(I,p,r)
#     Krn = createCommMat(rand(r,n));
#     Ay1 = (kronProd(Krn,Ip)); Ay2 = (kronProd(In,vec(B)));
#     Ay = kronProd(Iq, Ay1*Ay2);
#     Bx1 = (kronProd(Iq,Krn)); Bx2 = (kronProd(vec(A),Ir));
#     Bx = kronProd(Bx1*Bx2, Ip);
#     out = Ay*dA + Bx*dB
#     return out
# end
#
# function fdJ_kronABold(A,B, dA, dB)
#     n,q = size(A); p,r = size(B);
#     Iq = 1.0I(q); Ip = 1.0I(p); In = 1.0I(n); Ir = 1.0I(r);
#     Inq = Matrix{Float64}(I,n,q); Ipr = Matrix{Float64}(I,p,r)
#     Krn = createCommMat(rand(r,n));
#     Ay1 = kron(Krn,Ip); Ay2 = (kron(In,vec(B)));
#     # Ay_term2 = Ay1*Ay2;
#     Ay = kron(Iq, Ay1*Ay2);
#
#     Bx1 = kron(Iq,Krn); Bx2 = (kron(vec(A),Ir));
#     # Bx_term1 = Bx1*Bx2;
#     Bx = kron(Bx1*Bx2, Ip);
#     out = Ay*dA + Bx*dB
#     return out
# end
#
#
# dy1 = fdJ_kronAB2((tempA),tempB,dAtemp,dBtemp);
# dy2 = fdJ_kronABold((tempA),tempB,dAtemp,dBtemp);
# println("err = ", norm(dy1 - dy2))
# @btime fdJ_kronAB2((tempA),tempB,dAtemp,dBtemp);
# @btime fdJ_kronABold((tempA),tempB,dAtemp,dBtemp);
# ##
# n,q = size(tempA); p,r = size(tempB);
# Iq = 1.0I(q); Ip = 1.0I(p); In = 1.0I(n); Ir = 1.0I(r);
# Inq = Matrix{Float64}(I,n,q); Ipr = Matrix{Float64}(I,p,r);
# Krn = createCommMat(rand(r,n));
# Ay1 = collect(kronProd(Krn,Ip)); Ay2 = collect(kronProd(In,vec(tempB)[:,:]));
# # Ay1*Ay2
# # Ay_term2 = zeros(size(Ay1,1), size(Ay2,2));
# # @btime mul!(Ay_term2, Ay1, Ay2);
# Ay = kronProd(Iq, Ay1*Ay2);
#
# Bx1 = kronProd(Iq,Krn); Bx2 = kronProd(vec(tempA),Ir);
# Bx_term1 = zeros(size(Bx1,1), size(Bx2,2));
# mul!(Bx_term1,  )
# Bx = kronProd((kronProd(Iq,Krn)*kronProd(vec(A),Ir)), Ip);
# out = Ay*dA + Bx*dB
##
# convert(::Type{Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{var"#193#226"{Array{Float64,2},Array{Joint,1}},ForwardDiff.Dual{ForwardDiff.Tag{var"#201#258"{Array{Float64,2},var"#Fudiff1#225"{Array{Joint,1}}},Float64},Float64,10}},ForwardDiff.Dual{ForwardDiff.Tag{var"#201#258"{Array{Float64,2},var"#Fudiff1#225"{Array{Joint,1}}},Float64},Float64,10},10}}}, ::ForwardDiff.Dual{ForwardDiff.Tag{var"#201#258"{Array{Float64,2},var"#Fudiff1#225"{Array{Joint,1}}},Float64},Float64,10})


## sparse vs Dense
using SparseArrays

sz = 100; density = 0.1;
tA = sprand(Float64, sz, sz, density);
tB = sprand(Float64, sz, sz, density);

tAd = Array(tA); tBd = Array(tB);

@btime tA*tB;
@btime tAd*tBd;
