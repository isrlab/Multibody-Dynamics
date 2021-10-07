include("simulate.jl")

using Ipopt
using LazyArrays
using GenericLinearAlgebra
using SparseArrays

## trim function
function trim_kronLazy(j::Vector{Joint},GravityInInertial::Vector{Float64};
  ix::Vector{Integer}=zeros(Integer,20*(nB-1)+1),
  iy::Vector{Integer}=zeros(Integer,7*(nB-1)))

  nB = length(j) + 1; # number of bodies in the joint tree
  # ix: 0 indicates that the corresponding variable is free during optimization

  # iy: 0 indicates that the corresponding second derivative term needs to be zero
  # (specify indices for iy starting from second body)
  # iy = [1;4;5] means that the ẍ. β̈1 and β̈2 are free for j[1].RB2

  # Find which variables are free for optimization
  freeVarInd = findall(x->x==0,ix); # indices of free variables

  # find which acceleration-level terms are desired to be 0/constrained
  freeAccInd = findall(x->x==1,iy) # indices of free accelerations

  # Initial guesses
  # Initial guess is provided initial location
  x0 = []; # Initial guess
  for k = 1:length(j)
    append!(x0,j[k].RB2.x)
  end
  extFListCopy = externalFTotal(0.0,j);
  u0 = genU(extFListCopy); # Initial guess for control input

  n_xu = length(freeVarInd)-1; # Number of optimization variables
  lb = [-Inf*ones(n_xu);0.0]; ub = [Inf*ones(n_xu);Inf] # Bounds on optimization variables
  # lb[4] = 0.0;
  # lb[14+4] = 0.0;

  m = 1; # Number of nonlinear constraints (xdot = f(x,u)) for second order variables
  lb_g = [-Inf]; ub_g = [0.0];

  # extract states and control from optimization variables
  function getXU(z::AbstractArray{T}) where T<:Real
    z0 = [x0;u0[:,2:end][:];0.0];
    z0[freeVarInd] .= 0;

    Mz = zeros(length(z0), length(freeVarInd))
    for (i,ind) in enumerate(freeVarInd)
      Mz[ind,i] = 1;
    end

    zFull = z0 + Mz*z;

    X = zFull[1:14*(nB-1)];
    U = zFull[14*(nB-1)+1:end-1]
    return X,U
  end

  ## Constraint functions
  function qddFn(xuVar::AbstractArray{T}) where T<:Real
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    xdot = fxdot(x,u,j,g);
    y = Vector{T}(undef,7*nB); # qdd
    for i=1:nB
      y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
    end
    deleteat!(y,1:7); # inertial frame has no acceleration
    deleteat!(y,freeAccInd); # free acceleration variables
    return sparse(y)
  end

  function quatNormConstr(xuVar::AbstractArray{T}) where T<:Real
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    q_constr = Vector{T}(undef, nB)
    for i=1:nB
        q_constr[i] = norm(x[14*(i-1)+4:14*(i-1)+7])^2 - 1;
    end
    return q_constr
  end

  function jointLocConstr(xuVar::AbstractArray{T}) where T<:Real
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    jLoc_constr = Matrix{T}(undef,3,length(j))
    for i=1:length(j)
        b1_id = j[i].RB1.bodyID; b2_id = j[i].RB2.bodyID;
        x1 = x[(b1_id-1)*14+1:b1_id*14]
        x2 = x[(b2_id-1)*14+1:b2_id*14]
        j_b1 = x1[1:3] + transpose(quat2dcm(x1[4:7]))*j[i].pos1;
        j_b2 = x2[1:3] + transpose(quat2dcm(x2[4:7]))*j[i].pos2;
        jLoc_constr[:,i] = j_b1 - j_b2;
    end
    return jLoc_constr
  end

  function _nonlinear_constraints(z::AbstractArray{T}) where T<:Real
    xVar, uVar = getXU(z)
    x = [Vector{T}(j[1].RB1.x); xVar]
    u = [zeros(T,6) reshape(uVar,(6, nB-1))]
    xuVar = [x;u[:]]; # Vector of x,u

    # Dynamics constraint
    y = qddFn(xuVar);
    # Quaternion norm constraint
    q_constr = quatNormConstr(xuVar);
    # JointLocation constraint  (translation constraint)
    jLoc_constr = jointLocConstr(xuVar);

    constr = ([y;q_constr;jLoc_constr[:]]);

    J = sum(constr.*constr) - z[end];
    return J
  end

  ## Helper functions for derivatives
  function kronSum(A,B)
    # kronecker sum of 2 arrays
    out = kronProd(A, Matrix(I,(size(B)))) + kronProd(Matrix(I,(size(A))),B)
    return (out)
  end

  function kronProd(A,B)
    K = kron((A),(B));
    return sparse(K)
  end

  function createCommMat(M)
    # commutation matrix
    r,m = size(M);
    return createCommMat(r,m);
  end

  function createCommMat(r,m)
    # commutation matrix for a matrix of size r,m
    K = spzeros(m*r,m*r);
    Ir = sparse(I,r,r); Im = sparse(I,m,m);
    for i=1:r
      for j=1:m
          ei = @view Ir[:,i]; # ith-canonical unit vector of dimension r
          ej = @view Im[:,j];# jth-canonical unit vector of dimension m
          K += kronProd(ei*permutedims(ej), ej*permutedims(ei));
      end
    end
    return (K)
  end

  function fdJ_kronAB(A,B, dA, dB)
    n,q = size(A[:,:]); p,r = size(B[:,:]);
    Iq = sparse(1.0I,q,q); Ip = sparse(1.0I,p,p); In = sparse(1.0I,n,n); Ir = sparse(1.0I,r,r);
    Inq = sparse(1.0I,n,q); Ipr = sparse(1.0I,p,r)
    Krn = sparse(createCommMat(r,n));
    Ay1 = sparse(kronProd(Krn,Ip)); Ay2 = sparse(kronProd(In,vec(B)));
    Ay = sparse(kronProd(Iq, Ay1*Ay2));
    Bx1 = sparse(kronProd(Iq,Krn)); Bx2 = sparse(kronProd(vec(A),Ir));
    Bx = sparse(kronProd(Bx1*Bx2, Ip));
    out = sparse(Ay*dA + Bx*dB)
    return out
  end

  function grad_hess_nonlinear_constraints(z, flag)
    ## First derivatives
    xVar, uVar = getXU(z)
    xv = [j[1].RB1.x; xVar]
    uv = [zeros(6) reshape(uVar,(6, nB-1))];
    xu_var = [xv;uv[:]];

    function Mdiff1(x)
      # First derivatives of mass matrix and powers (1/2, -1/2, -1)
      M = Array(assembleM(x,j));
      Ms = (real(sqrt(M)));
      M_inv = (inv(M))
      Ms_inv = (real(sqrt(M_inv)));

      dM = sparse(ForwardDiff.jacobian(z -> assembleM(z,j),x))
      dM_inv = sparse(-(kronProd(permutedims(M_inv),M_inv))*dM)
      dMs = sparse(inv(Array(kronSum(permutedims(Ms), Ms)))*dM)
      dMs_inv = sparse(-kronProd(permutedims(Ms_inv),Ms_inv)*(dMs));
      return sparse(M), sparse(Ms), sparse(Ms_inv), dM, dMs, dMs_inv, sparse(M_inv), dM_inv
    end
    M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff1(xv);

    function Abdiff1(x)
      A, b = Ab_VecOfMat(x,j);
      dA = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],x);
      db = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[2],x);
      return sparse(A), sparse(dA), sparse(b), sparse(db)
    end
    A, dA, b, db = Abdiff1(xv);

    function Gpdiff1(x)
      function t1diff1(x)
        # t1 = A'*(AM)^{-1}*A';
        sz_A1 = size(A,1)
        t1 = sparse(permutedims(A)*inv(Array(A*M_inv*permutedims(A))))

        dA_tr = sparse(createCommMat(A)*dA);

        Y = A*M_inv*permutedims(A); invY = inv(Array(Y));
        dY1 = sparse(kronProd(A*M_inv,1.0I(sz_A1))*dA);
        dY2 = sparse(kronProd(A,A)*dM_inv);
        dY3 = sparse(kronProd(1.0I(sz_A1), A*M_inv)*dA_tr);
        dY = dY1 + dY2 + dY3;
        dY_inv = sparse(-kronProd(permutedims(invY),invY)*dY);

        dt1_1 = sparse(kronProd(permutedims(invY),1.0I(size(A,2)))*dA_tr);
        dt1_2 = sparse(kronProd(1.0I(sz_A1),permutedims(A))*dY_inv;)
        dt1 = sparse(dt1_1 + dt1_2);
        return t1, dt1
      end
      G = Array(A*Ms_inv);
      Gp = sparse(permutedims(G)*((G*permutedims(G))\I(size(G,1))))

      t1, dt1 = t1diff1(x)
      dGp1 = sparse(kronProd(permutedims(t1),1.0I(size(M,1)))*dMs_inv)

      sz_t2 =size(t1,2)
      t2 = sparse(kronProd(1.0I(sz_t2), Ms_inv))
      dGp2 = sparse(t2*dt1);

      dGp = sparse(dGp1 + dGp2)
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

    Fu = convert(SparseVector{Float64}, sparse(Fu));
    dFu = convert(SparseMatrixCSC{Float64}, sparse(dFu));
    dFu_u = convert(SparseMatrixCSC{Float64}, sparse(dFu_u));

    function hFn(x,u)
      h = b - A*M_inv*Fu;
      return sparse(h)
    end

    function hdiff1(x,u)
      # has derivatives wrt both x,u
      h = hFn(x,u);

      # dh
      term1_1 = sparse(kronProd(permutedims(M_inv*Fu), 1.0I(size(A,1))));
      # println("typeof(term1_1) = ", typeof(term1_1))
      # println("typeof(dA) = ", typeof(dA))
      term1 = sparse(term1_1*dA);
      term2 = sparse(kronProd(permutedims(Fu), A)*dM_inv);
      term3 = sparse(A*M_inv*dFu);
      dh_x = sparse(db - (term1 + term2 + term3));

      # dh_u
      dh_u = sparse(-A*M_inv*dFu_u);

      return sparse(h), sparse(dh_x), sparse(dh_u)
    end
    h, dh_x, dh_u = hdiff1(xv,uv);

    function Fcdiff1(x,u)
      Fc = Ms*Gp*h

      dFc1 = kronProd(permutedims(Gp*h),I(size(M,1)))*dMs
      dFc2 = kronProd(permutedims(h),Ms)*dGp
      dFc3 = Ms*Gp*dh_x;

      dFc = sparse(dFc1 + dFc2 + dFc3);

      dFc_u = Ms*Gp*dh_u;
      return sparse(Fc), sparse(dFc), sparse(dFc_u)
    end
    Fc, dFc, dFc_u = Fcdiff1(xv,uv);

    function accdiff1(x,u)
      # derivative of qdd w.r.t x
      dacc1 = sparse(kronProd(transpose(Fc+Fu),1.0I(size(M_inv,1)))*dM_inv)
      dacc2 = sparse(M_inv*(dFc + dFu))
      dqdd_x = sparse(dacc1 + dacc2)

      return (dqdd_x)
    end
    dqdd_x = accdiff1(xv,uv);
    ## Constraint derivative functions for gradient

    function dqddFn(xuVar)
      dqdd_u = sparse(M_inv*(dFc_u + dFu_u));
      out = [dqdd_x dqdd_u];
      return sparse(out)
    end
    qdd = qddFn(xu_var);
    qdd = convert(SparseVector{Float64}, sparse(qdd));
    dqdd = dqddFn(xu_var);

    function dynConstr_jac(xuVar)
      out = (2*permutedims(qdd)*dqdd)[:];
      return out
    end

    function quatNormConstr_diff(xuVar)
      f(y) = sum(quatNormConstr(y).*quatNormConstr(y));
      df(x) = ForwardDiff.gradient(y->f(y), x);
      grad_f = df(xuVar);
      d2f = ForwardDiff.jacobian(y->df(y), xuVar);
      return grad_f, d2f
    end

    function jointLocConstr_diff(xuVar)
      f(y) = sum((jointLocConstr(y)[:]).*(jointLocConstr(y)[:]));
      df(x) = ForwardDiff.gradient(y->f(y), x);
      grad_f = df(xuVar);
      d2f = ForwardDiff.jacobian(y->df(y), xuVar);
      return grad_f, d2f
    end

    function grad_nonlinear_constraints()
      dqdd_constr = dynConstr_jac(xu_var); # Dyn
      dquatNorm_constr,_ = quatNormConstr_diff(xu_var); # Quat
      djointLoc_constr, _ = jointLocConstr_diff(xu_var); # JointLoc

      dconstr_z = dqdd_constr + dquatNorm_constr + djointLoc_constr;

      # take care of free optimization variables
      freeVarInd_grad = 14 .+ freeVarInd[1:end-1];
      grad_constr_freeInd = dconstr_z[freeVarInd_grad];
      out = [grad_constr_freeInd;-1.0];
      return out
    end

    if flag == 1
      ## First derivatives required anyway
      out = grad_nonlinear_constraints();
    else
      ##  Derivatives for hessian
      function Mdiff2(x)
        # dM_xx
        dM_fn(y) = ForwardDiff.jacobian(z->assembleM(z,j),y);
        dM_xx = sparse(ForwardDiff.jacobian(z->dM_fn(z),x));

        # dMs_xx
        function z1Fn(x)
          z1= kronSum(permutedims(Ms), Ms);
          z1_inv = sparse(inv(Array(z1)));
          dMs_t = sparse(createCommMat(Ms)*dMs);
          dz1_1 = fdJ_kronAB(permutedims(Ms), sparse(1.0I, size(Ms)...), dMs_t, spzeros(size(dMs_t)...));
          dz1_2 = fdJ_kronAB(sparse(1.0I, size(Ms)...), Ms, spzeros(size(dMs)...), dMs);
          dz1 = sparse(dz1_1 + dz1_2);
          dz1_inv = sparse(-(kronProd(permutedims(z1_inv),z1_inv))*dz1);
          return z1_inv, dz1_inv
        end
        z1_inv, dz1_inv = z1Fn(x);
        out1 = sparse(kronProd(permutedims(dM), 1.0I(size(z1_inv,1)))*dz1_inv);
        out2 = sparse(kronProd(1.0I(size(dM,2)), z1_inv)*dM_xx);
        dMs_xx = sparse(out1 + out2);

        # dMsInv_xx
        R = kronProd(permutedims(Ms_inv), Ms_inv);
        dMs_inv_t = sparse(createCommMat(Ms_inv)*dMs_inv);
        dR = fdJ_kronAB(permutedims(Ms_inv), Ms_inv, dMs_inv_t, dMs_inv);

        out1 = sparse(kronProd(permutedims(dMs),sparse(1.0I(size(R,1))))*dR);
        out2 = kronProd(sparse(1.0I(size(dMs,2))),R)*dMs_xx;
        dMsInv_xx = sparse(-(out1 + out2));

        # dMinv_xx
        Q = kronProd(permutedims(M_inv), M_inv);
        dM_inv_t = sparse(createCommMat(M_inv)*dM_inv);
        dQ = fdJ_kronAB(permutedims(M_inv), M_inv, dM_inv_t, dM_inv);
        out1 = sparse(kronProd(permutedims(dM), sparse(1.0I(size(Q,1))))*dQ);
        out2 = sparse(kronProd(sparse(1.0I(size(dM,2))), Q)*dM_xx);
        dMinv_xx = sparse(-(out1 + out2));

        return dM_xx, dMs_xx, dMsInv_xx, dMinv_xx
      end
      dM_xx, dMs_xx, dMsInv_xx, dMinv_xx = Mdiff2(xv);

      function Abdiff2(x)

        dAFn(y) = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[1],y);
        cfg1 = ForwardDiff.JacobianConfig(dAFn, x, ForwardDiff.Chunk{1}());
        dA_xx = ForwardDiff.jacobian(dAFn,x, cfg1);
        dbFn(y) = ForwardDiff.jacobian(z->Ab_VecOfMat(z,j)[2],y);
        cfg1 = ForwardDiff.JacobianConfig(dbFn, x, ForwardDiff.Chunk{1}());
        db_xx = ForwardDiff.jacobian(dbFn,x, cfg1);
        return sparse(dA_xx), sparse(db_xx);
      end
      dA_xx, db_xx = Abdiff2(xv);

      function Gdiff2(x,u)
        G = sparse(A*Ms_inv);
        #dG
        dG1 = sparse(kronProd(permutedims(Ms_inv),sparse(1.0I(size(A,1))))*dA);
        dG2 = sparse(kronProd(sparse(1.0I(size(Ms_inv,2))),A)*dMs_inv);
        dG = sparse(dG1 + dG2);

        # dG_xx
        Q1 = kronProd(permutedims(Ms_inv), sparse(1.0I(size(A,1))));
        dMs_inv_t = sparse(createCommMat(Ms_inv)*dMs_inv);
        dQ1 = fdJ_kronAB(permutedims(Ms_inv), sparse(1.0I(size(A,1))),dMs_inv_t,spzeros(size(A,1)^2,length(x)));

        Q2 = kronProd(sparse(1.0I(size(Ms_inv,2))), A);
        dQ2 = fdJ_kronAB(sparse(1.0I(size(Ms_inv,2))), A, spzeros(size(dMs_inv)...), dA);

        out1 = sparse(kronProd(permutedims(dA),sparse(1.0I(size(Q1,1))))*dQ1);
        out2 = sparse(kronProd(sparse(1.0I(size(dA,2))), Q1)*dA_xx);
        out3 = sparse(kronProd(permutedims(dMs_inv), sparse(1.0I(size(Q2,1))))*dQ2);
        out4 = sparse(kronProd(sparse(1.0I(size(dMs_inv,2))), Q2)*dMsInv_xx);

        dG_xx = sparse(out1 + out2 + out3 + out4);
        return G, dG, dG_xx
      end
      G, dG, dG_xx = Gdiff2(xv,uv);
      # println("typeof(dG_xx) = ", typeof(dG_xx))

      function Hdiff2(x,u)
        # H = G*G^T
        dG_t = sparse(createCommMat(G)*dG);

        H = sparse(G*permutedims(G));

        dH = sparse(kronProd(G,sparse(1.0I(size(G,1))))*dG + kronProd(sparse(1.0I(size(G,1))),G)*dG_t);

        P1 = kronProd(G,sparse(1.0I(size(G,1))));
        dP1 = fdJ_kronAB(G, sparse(1.0I(size(G,1))), dG, spzeros(size(G,1)^2,length(x)));

        P2 = kronProd(sparse(1.0I(size(G,1))), G);
        dP2 = fdJ_kronAB(sparse(1.0I(size(G,1))), G, spzeros(size(G,1)^2,length(x)), dG);

        dG2_t = sparse(kronProd(sparse(1.0I(size(dG,2))), createCommMat(G))*dG_xx);

        out1 = sparse(kronProd(permutedims(dG), sparse(1.0I(size(P1,1))))*dP1);
        out2 = sparse(kronProd(sparse(1.0I(size(dG,2))), P1)*dG_xx);
        out3 = sparse(kronProd(permutedims(dG_t), sparse(1.0I(size(P2,1))))*dP2);
        out4 = sparse(kronProd(sparse(1.0I(size(dG_t,2))), P2)*dG2_t);

        dH_xx = sparse(out1 + out2 + out3 + out4);
        return H, dH, dH_xx
      end
      H, dH, dH_xx = Hdiff2(xv,uv);
      # println("typeof(dH_xx) = ", typeof(dH_xx))

      function Gpdiff2(x,u)
        dG_t = sparse(createCommMat(G)*dG);
        dG2_t = sparse(kronProd(sparse(1.0I(size(dG,2))), createCommMat(G))*dG_xx);

        H_inv = sparse(inv(Array(H)));
        dH_inv = sparse(-kronProd(permutedims(H_inv), H_inv)*dH);
        dH_inv_t = sparse(createCommMat(H_inv)*dH_inv);

        W = kronProd(permutedims(H_inv), H_inv);
        dW = fdJ_kronAB(permutedims(H_inv), H_inv, dH_inv_t, dH_inv);
        h1 = sparse(kronProd(permutedims(dH),sparse(1.0I(size(W,1))))*dW);
        h2 = sparse(kronProd(sparse(1.0I(size(dH,2))),W)*dH_xx);
        dH_inv2 = sparse(-(h1 + h2));

        V1 = kronProd(permutedims(H_inv), sparse(1.0I(size(G,2))));
        dV1 = fdJ_kronAB(permutedims(H_inv), sparse(1.0I(size(G,2))), dH_inv_t, spzeros(size(G,2)^2, length(x)));

        V2 = kronProd(sparse(1.0I(size(H_inv,2))), permutedims(G));
        dV2 = fdJ_kronAB(sparse(1.0I(size(H_inv,2))),permutedims(G), zeros(size(dH_inv_t)), dG_t);

        out1 = sparse(kronProd(permutedims(dG_t), sparse(1.0I(size(V1,1))))*dV1);
        out2 = sparse(kronProd(sparse(1.0I(size(dG_t,2))), V1)*dG2_t);
        out3 = sparse(kronProd(permutedims(dH_inv), sparse(1.0I(size(V2,1))))*dV2);
        out4 = sparse(kronProd(sparse(1.0I(size(dH_inv,2))), V2)*dH_inv2);

        dGp_xx = sparse(out1 + out2 + out3 + out4);

        return dGp_xx;
      end
      dGp_xx = Gpdiff2(xv,uv);
      # println("typeof(dGp_xx) = ", typeof(dGp_xx))

      function Fudiff2(x,u)
        dFu_xx = ForwardDiff.jacobian(y->Fudiff1(y,u)[2],x);
        dFu_ux = ForwardDiff.jacobian(y->Fudiff1(x,y)[2],u);
        dFu_xu = ForwardDiff.jacobian(y->Fudiff1(y,u)[3],x);
        dFu_uu = ForwardDiff.jacobian(y->Fudiff1(x,u)[3],u);
        return dFu_xx, dFu_ux, dFu_xu, dFu_uu
      end
      dFu_xx, dFu_ux, dFu_xu, dFu_uu = Fudiff2(xv,uv);

      dFu_xx = convert(SparseMatrixCSC{Float64}, sparse(dFu_xx));
      dFu_ux = convert(SparseMatrixCSC{Float64}, sparse(dFu_ux));
      dFu_xu = convert(SparseMatrixCSC{Float64}, sparse(dFu_xu));
      dFu_uu = convert(SparseMatrixCSC{Float64}, sparse(dFu_uu));
      Fu = convert(SparseVector{Float64}, sparse(Fu));
      dFu = convert(SparseMatrixCSC{Float64}, sparse(dFu));
      dFu_u = convert(SparseMatrixCSC{Float64}, sparse(dFu_u));

      function hdiff2(x,u)
        # dh_xx
        # terms 1,2
        T1 = kronProd(permutedims(M_inv*Fu), sparse(1.0I(size(A,1))))
        dMinvFu = sparse(kronProd(permutedims(Fu), sparse(1.0I(size(M_inv,1))))*dM_inv + M_inv*dFu);
        dMinvFu_t = dMinvFu;

        dT1 = fdJ_kronAB(permutedims(M_inv*Fu), sparse(1.0I(size(A,1))), dMinvFu_t, spzeros(size(A,1)^2, length(x)));
        term1 = sparse(kronProd(permutedims(dA), sparse(1.0I(size(T1,1))))*dT1)
        term2 = sparse(kronProd(sparse(1.0I(size(dA,2))), T1)*dA_xx);

        # terms 3,4
        T2 = kronProd(permutedims(Fu), A);
        dT2 = fdJ_kronAB(permutedims(Fu), A, dFu, dA);
        term3 = sparse(kronProd(permutedims(dM_inv), sparse(1.0I(size(T2,1))))*dT2);
        term4 = sparse(kronProd(sparse(1.0I(size(dM_inv,2))), T2)*dMinv_xx);

        # terms 5,6
        T3 = sparse(A*M_inv);
        dT3 = sparse(kronProd(permutedims(M_inv), sparse(1.0I(size(A,1))))*dA + kronProd(sparse(1.0I(size(M_inv,2))), A)*dM_inv);
        term5 = sparse(kronProd(permutedims(dFu), sparse(1.0I(size(A,1))))*dT3);
        term6 = sparse(kronProd(sparse(1.0I(length(x))), T3)*dFu_xx);

        dh_xx = sparse(db_xx - (term1 + term2 + term3 + term4 + term5 + term6));

        # dh_ux
        E1 = kronProd(permutedims(M_inv*Fu), sparse(1.0I(size(A,1))));
        dE1_u = fdJ_kronAB(permutedims(M_inv*Fu), sparse(1.0I(size(A,1))), M_inv*dFu_u, spzeros(size(A,1)^2, length(u)));

        E2 = kronProd(permutedims(Fu), A);
        dE2_u = fdJ_kronAB(permutedims(Fu), A, dFu_u, spzeros(length(A), length(u)));

        term1 = sparse(kronProd(permutedims(dA),sparse(1.0I(size(E1,1))))*dE1_u);
        term2 = sparse(kronProd(permutedims(dM_inv), sparse(1.0I(size(E2,1))))*dE2_u);
        term3 = sparse(kronProd(sparse(1.0I(size(dFu,2))), A*M_inv)*dFu_ux);

        dh_ux = sparse(-(term1 + term2 + term3));

        # dh_xu
        out1 = sparse(kronProd(permutedims(M_inv*dFu_u), sparse(1.0I(size(A,1))))*dA);
        out2 = sparse(kronProd(permutedims(dFu_u), A)*dM_inv);
        out3 = sparse(kronProd(sparse(1.0I(length(u))), A*M_inv)*dFu_xu);
        dh_xu = sparse(-(out1 + out2 + out3));

        # dh_uu
        dh_uu = sparse(-kronProd(sparse(1.0I(length(u))), A*M_inv)*dFu_uu);

        return dh_xx, dh_ux, dh_xu, dh_uu
      end
      dh_xx, dh_ux, dh_xu, dh_uu = hdiff2(xv,uv);

      h, dh_x, dh_u = hdiff1(xv,uv);

      function Fcdiff2(x,u)
        # dFc_xx
        function jacT1(x)
          D1 = kronProd(permutedims(Gp*h), sparse(1.0I(size(Ms,1))));
          dGph = sparse(kronProd(permutedims(h), sparse(1.0I(size(Gp,1))))*dGp + Gp*dh_x);
          dGph_t = dGph;
          dD1 = fdJ_kronAB(permutedims(Gp*h), sparse(1.0I(size(Ms,1))), dGph_t, spzeros(size(dMs)...));

          dMs_t = sparse(createCommMat(Ms)*dMs);

          dT1_1 = sparse(kronProd(permutedims(dMs), sparse(1.0I(size(D1,1))))*dD1);
          dT1_2 = sparse(kronProd(sparse(1.0I(size(dMs,2))), D1)*dMs_xx);

          dT1 = sparse(dT1_1 + dT1_2);
          return dT1
        end

        function jacT2(x)
          D2 = kronProd(permutedims(h), Ms);
          dD2 = fdJ_kronAB(permutedims(h), Ms, dh_x, dMs);
          # dGp_xx = Gpdiff2(x,u);

          dT2_1 = sparse(kronProd(permutedims(dGp),sparse(1.0I(size(D2,1))))*dD2);
          dT2_2 = sparse(kronProd(sparse(1.0I(size(dGp,2))), D2)*dGp_xx);
          dT2 = sparse(dT2_1 + dT2_2);
          return dT2
        end

        function jacT3(x)
          dT3_1 = sparse(kronProd(permutedims(Gp*dh_x), sparse(1.0I(size(Ms,1))))*dMs);
          dT3_2 = sparse(kronProd(permutedims(dh_x), Ms)*dGp);
          dT3_3 = sparse(kronProd(sparse(1.0I(size(dh_x,2))), Ms*Gp)*dh_xx);
          dT3 = sparse(dT3_1 + dT3_2 + dT3_3);

          return dT3
        end
        dFc_xx = sparse(jacT1(x) + jacT2(x) + jacT3(x));

        # dFc_ux
        # # first term
        D1 = kronProd(permutedims(Gp*h), sparse(1.0I(size(Ms,1))));
        dGph_u = sparse(Gp*dh_u);
        dGph_u_t = dGph_u;
        dD1_u = fdJ_kronAB(permutedims(Gp*h), sparse(1.0I(size(Ms,1))), dGph_u_t, spzeros(size(dMs,1),length(u)));
        out1 = sparse(kronProd(permutedims(dMs), sparse(1.0I(size(D1,1))))*dD1_u);
        # # second term
        D2 = sparse(kronProd(permutedims(h), Ms)*dGp);
        dD2_u = fdJ_kronAB(permutedims(h), Ms, dh_u, spzeros(size(dMs,1), length(u)));
        out2 = sparse(kronProd(permutedims(dGp), sparse(1.0I(size(D2,1))))*dD2_u);
        # # third term
        out3 = sparse(kronProd(sparse(1.0I(size(dh_x,2))), Ms*Gp)*dh_ux);
        dFc_ux = sparse(out1 + out2 + out3);

        # dFc_xu
        out1 = sparse(kronProd(permutedims(Gp*dh_u), sparse(1.0I(size(Ms,1))))*dMs);
        out2 = sparse(kronProd(permutedims(dh_u), Ms)*dGp);
        out3 = sparse(kronProd(sparse(1.0I(length(u))), Ms*Gp)*dh_xu);
        dFc_xu = sparse(out1 + out2 + out3);


        # dFc_uu
        dFc_uu = sparse(kronProd(sparse(1.0I(size(dh_u,2))), Ms*Gp)*dh_uu);

        return dFc_xx, dFc_ux, dFc_xu, dFc_uu
      end
      dFc_xx, dFc_ux, dFc_xu, dFc_uu = Fcdiff2(xv,uv);

      function hessqdd(x,u)
        T1 = kronProd(permutedims(Fc + Fu), sparse(1.0I(size(M_inv,1))));
        dT1 = fdJ_kronAB(permutedims(Fc + Fu), sparse(1.0I(size(M_inv,1))), (dFc + dFu), spzeros(size(dM_inv)...));
        T2 = sparse(dFc + dFu);
        dT2 = sparse(dFc_xx + dFu_xx);

        out1 = sparse(kronProd(permutedims(dM_inv), sparse(1.0I(size(T1,1))))*dT1);
        out2 = sparse(kronProd(sparse(1.0I(size(dM_inv,2))), T1)*dMinv_xx);
        out3 = sparse(kronProd(permutedims(T2), sparse(1.0I(size(M_inv,1))))*dM_inv);
        out4 = sparse(kronProd(sparse(1.0I(size(T2,2))), M_inv)*dT2);

        dqdd_xx = sparse(out1 + out2 + out3 + out4);
        return dqdd_xx
      end
      dqdd_xx = hessqdd(xv,uv);
      # println("typeof(dqdd_xx) = ", typeof(dqdd_xx))

      function cross_qdd_ux(x,u)
        T1 = kronProd(permutedims(Fc + Fu), sparse(1.0I(size(M_inv,1))));
        dT1_u = fdJ_kronAB(permutedims(Fc + Fu), sparse(1.0I(size(M_inv,1))), dFc_u + dFu_u, spzeros(size(M_inv,1)^2, length(u)));

        T2 = sparse(dFc + dFu);
        dT2_u = sparse(dFc_ux + dFu_ux);

        out1 = sparse(kronProd(permutedims(dM_inv), sparse(1.0I(size(T1,1))))*dT1_u);
        out2 = sparse(kronProd(sparse(1.0I(length(x))), M_inv)*dT2_u);
        dqdd_ux = sparse(out1+out2);
        return dqdd_ux
      end
      dqdd_ux = cross_qdd_ux(xv,uv);
      # println("typeof(dqdd_ux) = ", typeof(dqdd_ux))

      function cross_qdd_xu(x,u)
        out1 = sparse(kronProd(permutedims(dFc_u + dFu_u), sparse(1.0I(size(M_inv,1))))*dM_inv);
        out2 = sparse(kronProd(sparse(1.0I(length(u))), M_inv)*(dFc_xu + dFu_xu));
        dqdd_xu = sparse(out1 + out2);
        return dqdd_xu
      end
      dqdd_xu = cross_qdd_xu(xv,uv);
      # println("typeof(dqdd_xu) = ", typeof(dqdd_xu))

      function hessqdd_u(x,u)
        dqdd_uu = sparse(kronProd(sparse(1.0I(length(u))), M_inv)*(dFc_uu + dFu_uu));
        return dqdd_uu
      end
      dqdd_uu = hessqdd_u(xv,uv);
      # println("typeof(dqdd_uu) = ", typeof(dqdd_uu))

      function hessDyn(xuVar)
        x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end], (6,nB));

        y = qdd;
        dy = dqdd;

        d2y_xx = dqdd_xx;#  hessqdd(x,u);
        d2y_xu = dqdd_xu;# cross_qdd_xu(x,u);
        d2y_ux = dqdd_ux;# cross_qdd_ux(x,u);
        d2y_uu = dqdd_uu;# hessqdd_u(x,u);

        d2y = [d2y_xx d2y_ux; d2y_xu d2y_uu]
        # d2y = convert(SparseMatrixCSC{Float64},sparse([d2y_xx d2y_ux; d2y_xu d2y_uu]));

        out1 = sparse(kronProd(sparse(1.0I(length(xuVar))), permutedims(y))*d2y);
        # out1 = convert(SparseMatrixCSC{Float64}, out1);
        out2 = sparse(permutedims(dy)*dy);
        out = sparse(2*(out1 + out2)); # hessian of dynamics cost = sum(y.*y)
        return out;
      end
      dqdd_xuVar = hessDyn(xu_var);
      # println("typeof(dqdd_xuVar) = ", typeof(dqdd_xuVar))

      function hess_nonlinear_constraints()
        zLen = length(z);

        hess_dyn = hessDyn(xu_var); # hessian of dynamics cost
        _, hess_qNorm = quatNormConstr_diff(xu_var); # hessian of quatNormConstr cost
        _, hess_jLoc = jointLocConstr_diff(xu_var); # hessian of jointLocConstr cost

        hess_constr = hess_dyn +  hess_qNorm + hess_jLoc;

        # taking only free variables into consideration
        freeVarInd_hess = 14 .+ freeVarInd[1:end-1];
        hess_constr_free = hess_constr[freeVarInd_hess, freeVarInd_hess];

        out = zeros(typeof(hess_constr_free[1]),(zLen, zLen));
        out[1:zLen-1, 1:zLen - 1] = hess_constr_free;
        return out
      end
      out = hess_nonlinear_constraints();
    end
    return out
  end

  ## IpOpt functions
  function eval_f(z)
    J = z[end]
    return J
  end

  function eval_g(z, g)
    g[:] .= _nonlinear_constraints(z)
  end

  function eval_grad_f(z, grad_f)
    v = zeros(length(z)); v[end] = 1.0;
    grad_f[:] = v;
  end


  function eval_jac_g(z, mode, rows, cols, values)
    nvar = length(z);
    ng = 1;

    if mode == :Structure
      index = 1;
      for i=1:ng
          for j=1:nvar
              rows[index] = i;
              cols[index] = j;
              index += 1;
          end
      end
    else
      flag = 1;
      values[:] = grad_hess_nonlinear_constraints(z,flag);
    end
  end

  function eval_h(z, mode, rows, cols, obj_factor, lambda, values) # Hessian of cost.
    nvar = length(z);
    if mode == :Structure
      # Symmetric matrix, fill the lower left triangle only
      idx = 1
      for row = 1:nvar
        for col = 1:row
          rows[idx] = row
          cols[idx] = col
          idx += 1
        end
      end
    else
      # H = obj_factor*zeros(nvar, nvar)
      # for i=1:length(lambda)
      flag = 2;
      H = lambda[1]*grad_hess_nonlinear_constraints(z,flag);
      # end

      # Copy them into return variable -- only lower triangle.
      idx = 1;
      for row = 1:nvar
        for col = 1:row
          values[idx] = H[row,col];
          idx += 1;
        end
      end
    end
  end

  # Initial guess from freeVarInd
  n = length(freeVarInd);
  guessFull = [x0;u0[:,2:end][:];0.0];
  inGuess = guessFull[freeVarInd];

  nele_jac = Integer(m*n);
  nele_hess = Integer(n*(n+1)/2);
  prob = createProblem(n, lb, ub, m, lb_g, ub_g, nele_jac, nele_hess, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h);
  # Ipopt.addOption(prob, "linear_solver", "pardiso");
  prob.x = inGuess; # Initial guess.
  status = solveProblem(prob);

  trim_x = prob.x[1:length(x0)];
  gam = prob.x[end];
  trim_u = prob.x[length(x0)+1:end-1];
  return trim_x, trim_u, gam
end

## setup for pendulum
# m = 1.0; l = 1.0; # Mass and length of bar
# # Assuming bar revolves about Y axis
# I1 = [1 0 0; 0 m * l^2 / 12 0; 0 0 1];
# # Testing Dynamics with Revolute Joint
# R1 = RigidBody(m, I1, 2);
# RbI = InertialFrameAsRB()
#
# # Suspended at an angle theta from vertical
# theta = deg2rad(30)
# x_temp = [l/2*sin(theta); 0.0; -l/2*cos(theta)];
# x0R1 = ([x_temp;[1;zeros(3)];zeros(3);zeros(4)]);
# initialiseRigidBody!(R1,x0R1)
#
# # Axis about which bar is revolving
# axisY = [0.0 1.0 0.0][:];
#
# rj1 = [0.0 0.0 0.0][:]; # Joint Location in body frame of first body
# rj2 = -R1.x[1:3];
#
# j1 = Joint(RbI, R1, rj1, rj2, type="Revolute", axis=axisY);
# j = [j1]; # Joint tree for pendulum
# g = [0.0,0.0,-9.806]; # Gravity Vector.
# x0Orig, u0Orig = getXU_0(j); # Initial values
# nB = length(j) + 1; # number of bodies
# ix= zeros(Integer,20*(nB-1)+1); # indicates which variables are free for optimization
# ix[14*(nB-1)+1:20*(nB-1)] .= 1; ## keeping u constant, i.e., no force applied on pendulum
# trim_U = u0Orig;
# iy = zeros(Integer,7*(nB-1)); ## indicates which acceleration level terms need not be constrained. iy == 0 means acceleration is zero for corresponding state.
# out = trim_kronLazy(j,g,ix=ix, iy=iy);
# trim_x, trim_u, gam = out;
# println("trim_x = ", trim_x)
# trim_X = [j[1].RB1.x;trim_x]
# println("qconstr = ", norm(trim_x[4:7])-1)
# println("ẍ =", norm(fxdot(trim_X,trim_U, j, g)))

## setup for quadrotor
include("nRotor.jl");
j = gen_nRotor(4);
g =([0.0,0.0,-9.806]);
x0Orig, u0Orig = getXU_0(j); # Initial values
nB = length(j) + 1; # number of bodies
ix= zeros(Integer,20*(nB-1)+1); # indicates which variables are free for optimization
ix[14*(nB-1)+1:20*(nB-1)] .= 1; ## keeping u constant, i.e., hover
trim_U = u0Orig;
iy = zeros(Integer,7*(nB-1)); ## indicates which acceleration level terms need not be constrained. iy == 0 means acceleration is zero for corresponding state.
out = trim_kronLazy(j,g,ix=ix, iy=iy);
trim_x, trim_u, gam = out;
println("trim_x = ", trim_x)
trim_X = [j[1].RB1.x;trim_x]

println("qconstr = ", norm(trim_x[4:7])-1)
println("ẍ =", norm(fxdot(trim_X,trim_U, j, g)))
