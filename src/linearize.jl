include("simulate.jl")

using ForwardDiff
using FiniteDifferences
using LazyArrays

## helper functions for matrix derivatives
function kronSum(A,B)
  # kronecker sum of 2 arrays
  out = kronProd(A, Matrix(I,(size(B)))) + kronProd(Matrix(I,(size(A))),B)
  return out
end

function kronProd(A,B)
    K = ApplyArray(kron,A[:,:],B[:,:]);
    out = Matrix{Float64}(undef,size(K)...); copyto!(out,K);
    return out
end

function createCommMat(M)
  # commutation matrix
  r,m = size(M);
  K = zeros(m*r,m*r);
  for i=1:r
      for j=1:m
          ei = 1.0I(r)[:,i] # ith-canonical unit vector of dimension r
          ej = 1.0I(m)[:,j] # jth-canonical unit vector of dimension m
          K += kronProd(ei*permutedims(ej), ej*permutedims(ei))
      end
  end
  return K
end

## main function for linearize
function linearize(xv,uv,j,GravityInInertial)
    function Mdiff1(x)
        # First derivatives of mass matrix and powers (1/2, -1/2, -1)
        M = assembleM(x,j)
        Ms = real(sqrt(M))
        M_inv = inv(M)
        Ms_inv = real(sqrt(M_inv))

        dM = ForwardDiff.jacobian(z -> assembleM(z,j),x)
        dM_inv = -(kronProd(permutedims(M_inv),M_inv))*dM;
        dMs = inv(kronSum(permutedims(Ms), Ms))*dM
        dMs_inv = -kronProd(permutedims(Ms_inv),Ms_inv)*(dMs);
        return M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv
    end
    M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff1(xv);

    function Adiff1(x)
        A, _ = Ab_VecOfMat(x,j);
        dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
        return A, dA
    end
    A, dA = Adiff1(xv);

    function bdiff1(x)
        _, b = Ab_VecOfMat(x,j);
        db = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[2],x);
        return b, db
    end
    b, db = bdiff1(xv);

    function Gpdiff1(x)
        function t1diff1(x)
            # t1 = A'*(AM)^{-1}*A';
            sz_A1 = size(A,1)
            t1 = permutedims(A)*inv(A*inv(M)*permutedims(A))

            dA_tr = createCommMat(A)*dA;

            Y = A*M_inv*permutedims(A); invY = inv(Y);
            dY1 = kronProd(A*M_inv,1.0I(sz_A1))*dA;
            dY2 = kronProd(A,A)*dM_inv;
            dY3 = kronProd(1.0I(sz_A1), A*M_inv)*dA_tr;
            dY = dY1 + dY2 + dY3;
            dY_inv = -kronProd(permutedims(invY),invY)*dY;

            dt1_1 = kronProd(permutedims(invY),1.0I(size(A,2)))*dA_tr;
            dt1_2 = kronProd(1.0I(sz_A1),permutedims(A))*dY_inv;
            dt1 = dt1_1 + dt1_2;
            return t1, dt1
        end
        G = A*Ms_inv;
        Gp = permutedims(G)*((G*permutedims(G))\I(size(G,1)))

        t1, dt1 = t1diff1(x)
        dGp1 = kronProd(permutedims(t1),1.0I(size(M,1)))*dMs_inv

        sz_t2 =size(t1,2)
        t2 = kronProd(1.0I(sz_t2), Ms_inv)
        dGp2 = t2*dt1;

        dGp = dGp1 + dGp2
    return Gp, dGp
    end
    Gp, dGp = Gpdiff1(xv);

    function Fudiff1(x,u)
        Fu = assembleF(x,u,j,g);
        dFu = ForwardDiff.jacobian(y->assembleF(y,u,j,g),x);
        dFu_u = ForwardDiff.jacobian(y->assembleF(x,y,j,g),u);
        return Fu, dFu, dFu_u
    end
    Fu, dFu, dFu_u = Fudiff1(xv,uv);

    function hdiff1(x,u)
        # has derivatives wrt both x,u
        h = b - A*inv(M)*Fu;

        # dh
        term1_1 = kronProd(permutedims(M_inv*Fu), 1.0I(size(A,1)));
        term1 = term1_1*dA;
        term2 = kronProd(permutedims(Fu), A)*dM_inv;
        term3 = A*M_inv*dFu;
        dh_x = db - (term1 + term2 + term3);

        # dh_u
        dh_u = -A*M_inv*dFu_u;

        return h, dh_x, dh_u
    end
    h, dh_x, dh_u = hdiff1(xv,uv);

    function Fcdiff1(x,u)
        Fc = Ms*Gp*h

        dFc1 = kronProd(permutedims(Gp*h),I(size(M,1)))*dMs
        dFc2 = kronProd(permutedims(h),Ms)*dGp
        dFc3 = Ms*Gp*dh_x;

        dFc = dFc1 + dFc2 + dFc3;

        dFc_u = Ms*Gp*dh_u;
        return Fc, dFc, dFc_u
    end
    Fc, dFc, dFc_u = Fcdiff1(xv,uv);

    function accdiff1(x,u)
        # derivative of qdd w.r.t x
        dacc1 = kronProd(permutedims(Fc+Fu),1.0I(size(M_inv,1)))*dM_inv
        dacc2 = M_inv*(dFc + dFu)
        dqdd_x = dacc1 + dacc2

        return dqdd_x
    end
    dqdd_x = accdiff1(xv,uv);
    dqdd_u = M_inv*(dFc_u + dFu_u);

    function fxdotDiff(x,u)
        len_x = length(x); len_u = length(u);
        df_x = zeros(len_x,len_x)

        df_u = zeros(len_x, len_u);

        for k=1:length(j)
            df_x[14*k+1:14*k+7,14*k+8:14*k+14] = 1.0*I(7)
            df_x[14*k+8:14*k+14,:] = dqdd_x[7*(k-1)+1:7*k,:]

            df_u[14*k+8:14*k+14,:] = dqdd_u[7*(k-1)+1:7*k,:];
        end
        df_x_woIn = df_x[15:end,15:end]; # without the inertial frame
        df_u_woIn = df_u[15:end,7:end];  # without the inertial frame
        return df_x_woIn, df_u_woIn
    end

    dfx_x, dfx_u = fxdotDiff(xv,uv);
    return dfx_x, dfx_u
end
## Test for different robots
# x0Orig, u0Orig = getXU_0(j);
# fdJ_A, fdJ_B = linearize(x0Orig, u0Orig, j, g);
#
# m = central_fdm(5,1); # first order differential using fifth order central difference method
# finJ_A = FiniteDifferences.jacobian(m, z->fxdot(z,u0Orig,j,g),x0Orig)[1];
# finJ_B = FiniteDifferences.jacobian(m, z->fxdot(x0Orig,z,j,g),u0Orig)[1];
#
# println("err_A = ", norm(fdJ_A - finJ_A))
# println("err_B = ", norm(fdJ_B - finJ_B))
