# include("simulate.jl")

using Ipopt
using GenericLinearAlgebra

## Ipopt directly
function trim_kron(j::Vector{Joint},GravityInInertial::Vector{Float64};
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
  function qddFn(xuVar)
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    xdot = fxdot(x,u,j,g);
    y = Vector(undef,7*nB); # qdd
    for i=1:nB
      y[7*(i-1)+1:7*(i)] = xdot[7*(2*i-1)+1:14*i]
    end
    deleteat!(y,1:7); # inertial frame has no acceleration
    deleteat!(y,freeAccInd); # free acceleration variables
    return y
  end

  function quatNormConstr(xuVar)
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    q_constr = Vector(undef, nB)
    for i=1:nB
        q_constr[i] = norm(x[14*(i-1)+4:14*(i-1)+7])^2 - 1;
    end
    return q_constr
  end

  function jointLocConstr(xuVar)
    x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end],(6,nB));
    jLoc_constr = Matrix(undef,3,length(j))
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

    constr = [y;q_constr;jLoc_constr[:]];

    J = sum(constr.*constr) - z[end];
    return J
  end

  ## Helper functions for derivatives
  function kronSum(A,B)
    # kronecker sum of 2 arrays
    out = kron(A, Matrix(I,(size(B)))) + kron(Matrix(I,(size(A))),B)
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
            K += kron(ei*permutedims(ej), ej*permutedims(ei))
        end
    end
    return K
  end

  function fdJ_kronAB(A,B, dA, dB)
      n,q = size(A); p,r = size(B);
      Iq = 1.0I(q); Ip = 1.0I(p); In = 1.0I(n); Ir = 1.0I(r);
      Inq = Matrix{Float64}(I,n,q); Ipr = Matrix{Float64}(I,p,r)
      Krn = createCommMat(rand(r,n));
      Ay1 = (kron(Krn,Ip)); Ay2 = (kron(In,vec(B)));
      Ay = kron(Iq, Ay1*Ay2);
      Bx = kron((kron(Iq,Krn)*kron(vec(A),Ir)), Ip);
      out = Ay*dA + Bx*dB
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
      M = assembleM(x,j)
      Ms = real(sqrt(M))
      M_inv = inv(M)
      Ms_inv = real(sqrt(M_inv))

      dM = ForwardDiff.jacobian(z -> assembleM(z,j),x)
      dM_inv = -(kron(permutedims(M_inv),M_inv))*dM;
      dMs = inv(kronSum(permutedims(Ms), Ms))*dM
      dMs_inv = -kron(permutedims(Ms_inv),Ms_inv)*(dMs);
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
      G = A*Ms_inv;
      Gp = permutedims(G)*((G*permutedims(G))\I(size(G,1)))

      t1, dt1 = t1diff1(x)
      dGp1 = kron(permutedims(t1),1.0I(size(M,1)))*dMs_inv

      sz_t2 =size(t1,2)
      t2 = kron(1.0I(sz_t2), Ms_inv)
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

    function hFn(x,u)
      A,b = Ab_VecOfMat(x,j);
      Fu = assembleF(x,u,j,g);
      M = assembleM(x,j);
      h = b - A*inv(M)*Fu;
      return h
    end

    function hdiff1(x,u)
      # has derivatives wrt both x,u
      Fu, dFu, dFu_u = Fudiff1(x,u);
      A, dA = Adiff1(x);
      b, db = bdiff1(x);
      M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff1(x);

      h = b - A*inv(M)*Fu;

      # dh
      term1 = kron(permutedims(M_inv*Fu), 1.0I(size(A,1)))*dA;
      term2 = kron(permutedims(Fu), A)*dM_inv;
      term3 = A*M_inv*dFu;
      dh_x = db - (term1 + term2 + term3);

      # dh_u
      dh_u = -A*M_inv*dFu_u;

      return h, dh_x, dh_u
    end
    h, dh_x, dh_u = hdiff1(xv,uv);

    function Fcdiff1(x,u)
      Fc = Ms*Gp*h

      dFc1 = kron(permutedims(Gp*h),I(size(M,1)))*dMs
      dFc2 = kron(permutedims(h),Ms)*dGp
      dFc3 = Ms*Gp*dh_x;

      dFc = dFc1 + dFc2 + dFc3;

      dFc_u = Ms*Gp*dh_u;
      return Fc, dFc, dFc_u
    end
    Fc, dFc, dFc_u = Fcdiff1(xv,uv);

    function accdiff1(x,u)
      # derivative of qdd w.r.t x
      dacc1 = kron(transpose(Fc+Fu),1.0I(size(M_inv,1)))*dM_inv
      dacc2 = M_inv*(dFc + dFu)
      dqdd_x = dacc1 + dacc2

      return dqdd_x
    end
    dqdd_x = accdiff1(xv,uv);
    ## Constraint derivative functions for gradient

    function dqddFn(xuVar)
      dqdd_u = M_inv*(dFc_u + dFu_u);
      out = [dqdd_x dqdd_u];
      return out
    end
    qdd = qddFn(xu_var);
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
        M, Ms, Ms_inv, dM, dMs, dMs_inv, M_inv, dM_inv = Mdiff1(x);

        # dM_xx
        dM_fn(y) = ForwardDiff.jacobian(z->assembleM(z,j),y);
        dM_xx = ForwardDiff.jacobian(z->dM_fn(z),x);

        # dMs_xx
        function z1Fn(x)
          z1= kronSum(permutedims(Ms), Ms);
          z1_inv = inv(z1);
          dMs_t = createCommMat(Ms)*dMs;
          dz1_1 = fdJ_kronAB(permutedims(Ms),1.0I(size(Ms,1)),dMs_t,zeros(size(dMs_t)));
          dz1_2 = fdJ_kronAB(1.0I(size(Ms,2)),Ms,zeros(size(dMs)),dMs);
          dz1 = dz1_1 + dz1_2;
          dz1_inv = -(kron(permutedims(z1_inv),z1_inv))*dz1;
          return z1_inv, dz1_inv
        end
        z1_inv, dz1_inv = z1Fn(x);
        out1 = kron(permutedims(dM), 1.0I(size(z1_inv,1)))*dz1_inv;
        out2 = kron(1.0I(size(dM,2)), z1_inv)*dM_xx
        dMs_xx = out1 + out2;

        # dMsInv_xx
        R = kron(permutedims(Ms_inv), Ms_inv);
        dMs_inv_t = createCommMat(Ms_inv)*dMs_inv;
        dR = fdJ_kronAB(permutedims(Ms_inv), Ms_inv, dMs_inv_t, dMs_inv);

        out1 = kron(permutedims(dMs),1.0I(size(R,1)))*dR;
        out2 = kron(1.0I(size(dMs,2)),R)*dMs_xx;
        dMsInv_xx = -(out1 + out2);

        # dMinv_xx
        Q = kron(permutedims(M_inv), M_inv);
        dM_inv_t = createCommMat(M_inv)*dM_inv;
        dQ = fdJ_kronAB(permutedims(M_inv), M_inv, dM_inv_t, dM_inv);
        out1 = kron(permutedims(dM), 1.0I(size(Q,1)))*dQ;
        out2 = kron(1.0I(size(dM,2)), Q)*dM_xx;
        dMinv_xx = -(out1 + out2);

        return dM_xx, dMs_xx, dMsInv_xx, dMinv_xx
      end
      dM_xx, dMs_xx, dMsInv_xx, dMinv_xx = Mdiff2(xv);

      function Adiff2(x)
        dAFn(y) = Adiff1(y)[2];
        cfg1 = ForwardDiff.JacobianConfig(dAFn, x, ForwardDiff.Chunk{1}());
        dA_xx = ForwardDiff.jacobian(dAFn,x, cfg1);
        return dA_xx;
      end
      dA_xx = Adiff2(xv);
      A, dA = Adiff1(xv);

      function bdiff2(x)
        dbFn(y) = bdiff1(y)[2];
        cfg1 = ForwardDiff.JacobianConfig(dbFn, x, ForwardDiff.Chunk{1}());
        db_xx = ForwardDiff.jacobian(dbFn,x, cfg1);
        return db_xx;
      end
      db_xx = bdiff2(xv);
      b, db = bdiff1(xv);

      function Gdiff2(x,u)
        G = A*Ms_inv;
        #dG
        dG1 = kron(permutedims(Ms_inv),1.0I(size(A,1)))*dA;
        dG2 = kron(1.0I(size(Ms_inv,2)),A)*dMs_inv;
        dG = dG1 + dG2;

        # dG_xx
        Q1 = kron(permutedims(Ms_inv), 1.0I(size(A,1)));
        dMs_inv_t = createCommMat(Ms_inv)*dMs_inv;
        dQ1 = fdJ_kronAB(permutedims(Ms_inv), 1.0I(size(A,1)),dMs_inv_t,zeros(size(A,1)^2,length(x)));

        Q2 = kron(1.0I(size(Ms_inv,2)), A);
        dQ2 = fdJ_kronAB(1.0I(size(Ms_inv,2)), A, zeros(size(dMs_inv)), dA);

        out1 = kron(permutedims(dA),1.0I(size(Q1,1)))*dQ1;
        out2 = kron(1.0I(size(dA,2)), Q1)*dA_xx;
        out3 = kron(permutedims(dMs_inv), 1.0I(size(Q2,1)))*dQ2;
        out4 = kron(1.0I(size(dMs_inv,2)), Q2)*dMsInv_xx;

        dG_xx = out1 + out2 + out3 + out4;
        return G, dG, dG_xx
      end
      G, dG, dG_xx = Gdiff2(xv,uv);

      function Hdiff2(x,u)
        # H = G*G^T
        dG_t = createCommMat(G)*dG;

        H = G*permutedims(G);

        dH = kron(G,1.0I(size(G,1)))*dG + kron(1.0I(size(G,1)),G)*dG_t;

        P1 = kron(G,1.0I(size(G,1)));
        dP1 = fdJ_kronAB(G, 1.0I(size(G,1)), dG, zeros(size(G,1)^2,length(x)));

        P2 = kron(1.0I(size(G,1)), G);
        dP2 = fdJ_kronAB(1.0I(size(G,1)), G, zeros(size(G,1)^2,length(x)), dG);

        dG2_t = kron(1.0I(size(dG,2)), createCommMat(G))*dG_xx;

        out1 = kron(permutedims(dG), 1.0I(size(P1,1)))*dP1;
        out2 = kron(1.0I(size(dG,2)), P1)*dG_xx;
        out3 = kron(permutedims(dG_t), 1.0I(size(P2,1)))*dP2;
        out4 = kron(1.0I(size(dG_t,2)), P2)*dG2_t;

        dH_xx = out1 + out2 + out3 + out4;
        return H, dH, dH_xx
      end
      H, dH, dH_xx = Hdiff2(xv,uv);

      function Gpdiff2(x,u)
        dG_t = createCommMat(G)*dG;
        dG2_t = kron(1.0I(size(dG,2)), createCommMat(G))*dG_xx;

        H_inv = inv(H);
        dH_inv = -kron(permutedims(H_inv), H_inv)*dH;
        dH_inv_t = createCommMat(H_inv)*dH_inv;

        W = kron(permutedims(H_inv), H_inv);
        dW = fdJ_kronAB(permutedims(H_inv), H_inv, dH_inv_t, dH_inv);
        h1 = kron(permutedims(dH),1.0I(size(W,1)))*dW;
        h2 = kron(1.0I(size(dH,2)),W)*dH_xx;
        dH_inv2 = -(h1 + h2);

        V1 = kron(permutedims(H_inv), 1.0I(size(G,2)));
        dV1 = fdJ_kronAB(permutedims(H_inv), 1.0I(size(G,2)), dH_inv_t, zeros(size(G,2)^2, length(x)));

        V2 = kron(1.0I(size(H_inv,2)), permutedims(G));
        dV2 = fdJ_kronAB(1.0I(size(H_inv,2)),permutedims(G), zeros(size(dH_inv_t)), dG_t);

        out1 = kron(permutedims(dG_t), 1.0I(size(V1,1)))*dV1;
        out2 = kron(1.0I(size(dG_t,2)), V1)*dG2_t;
        out3 = kron(permutedims(dH_inv), 1.0I(size(V2,1)))*dV2;
        out4 = kron(1.0I(size(dH_inv,2)), V2)*dH_inv2;

        dGp_xx = out1 + out2 + out3 + out4;

        return dGp_xx;
      end
      dGp_xx = Gpdiff2(xv,uv);

      function Fudiff2(x,u)
        dFu_xx = ForwardDiff.jacobian(y->Fudiff1(y,u)[2],x);
        dFu_ux = ForwardDiff.jacobian(y->Fudiff1(x,y)[2],u);
        dFu_xu = ForwardDiff.jacobian(y->Fudiff1(y,u)[3],x);
        dFu_uu = ForwardDiff.jacobian(y->Fudiff1(x,u)[3],u);
        return dFu_xx, dFu_ux, dFu_xu, dFu_uu
      end
      dFu_xx, dFu_ux, dFu_xu, dFu_uu = Fudiff2(xv,uv);

      function hdiff2(x,u)
        # dh_xx
        # terms 1,2
        T1 = kron(permutedims(M_inv*Fu), 1.0I(size(A,1)))
        dMinvFu = kron(permutedims(Fu), 1.0I(size(M_inv,1)))*dM_inv + M_inv*dFu;
        dMinvFu_t = dMinvFu;
        dT1 = fdJ_kronAB(permutedims(M_inv*Fu), 1.0I(size(A,1)), dMinvFu_t, zeros(size(A,1)^2, length(x)));
        term1 = kron(permutedims(dA), 1.0I(size(T1,1)))*dT1
        term2 = kron(1.0I(size(dA,2)), T1)*dA_xx;

        # terms 3,4
        T2 = kron(permutedims(Fu), A);
        dT2 = fdJ_kronAB(permutedims(Fu), A, dFu, dA);
        term3 = kron(permutedims(dM_inv), 1.0I(size(T2,1)))*dT2;
        term4 = kron(1.0I(size(dM_inv,2)), T2)*dMinv_xx;

        # terms 5,6
        T3 = A*M_inv;
        dT3 = kron(permutedims(M_inv), 1.0I(size(A,1)))*dA + kron(1.0I(size(M_inv,2)), A)*dM_inv;
        term5 = kron(permutedims(dFu), 1.0I(size(A,1)))*dT3;
        term6 = kron(1.0I(length(x)), T3)*dFu_xx;

        dh_xx = db_xx - (term1 + term2 + term3 + term4 + term5 + term6);

        # dhFn(y) = ForwardDiff.jacobian(z->hFn(z,u),y);
        # cfg1 = ForwardDiff.JacobianConfig(dhFn, x, ForwardDiff.Chunk{1}());
        # dh_xx = ForwardDiff.jacobian(dhFn,x,cfg1);

        # dh_ux
        E1 = kron(permutedims(M_inv*Fu), 1.0I(size(A,1)));
        dE1_u = fdJ_kronAB(permutedims(M_inv*Fu), 1.0I(size(A,1)), M_inv*dFu_u, zeros(size(A,1)^2, length(u)));

        E2 = kron(permutedims(Fu), A);
        dE2_u = fdJ_kronAB(permutedims(Fu), A, dFu_u, zeros(length(A), length(u)));

        term1 = kron(permutedims(dA),1.0I(size(E1,1)))*dE1_u;
        term2 = kron(permutedims(dM_inv), 1.0I(size(E2,1)))*dE2_u;
        term3 = kron(1.0I(size(dFu,2)), A*M_inv)*dFu_ux;

        dh_ux = -(term1 + term2 + term3);

        # dh_xu
        out1 = kron(permutedims(M_inv*dFu_u), 1.0I(size(A,1)))*dA;
        out2 = kron(permutedims(dFu_u), A)*dM_inv;
        out3 = kron(1.0I(length(u)), A*M_inv)*dFu_xu;
        dh_xu = -(out1 + out2 + out3);

        # dh_uu
        dh_uu = -kron(1.0I(length(u)), A*M_inv)*dFu_uu;

        return dh_xx, dh_ux, dh_xu, dh_uu
      end
      dh_xx, dh_ux, dh_xu, dh_uu = hdiff2(xv,uv);
      h, dh_x, dh_u = hdiff1(xv,uv);

      function Fcdiff2(x,u)
        # dFc_xx
        function jacT1(x)
          D1 = kron(permutedims(Gp*h), 1.0I(size(Ms,1)));
          dGph = kron(permutedims(h), 1.0I(size(Gp,1)))*dGp + Gp*dh_x;
          dGph_t = dGph;
          dD1 = fdJ_kronAB(permutedims(Gp*h), 1.0I(size(Ms,1)), dGph_t, zeros(size(dMs)));

          dMs_t = createCommMat(Ms)*dMs;

          dT1_1 = kron(permutedims(dMs), 1.0I(size(D1,1)))*dD1;
          dT1_2 = kron(1.0I(size(dMs,2)), D1)*dMs_xx;

          dT1 = dT1_1 + dT1_2;
          return dT1
        end

        function jacT2(x)
          D2 = kron(permutedims(h), Ms);
          dD2 = fdJ_kronAB(permutedims(h), Ms, dh_x, dMs);
          # dGp_xx = Gpdiff2(x,u);

          dT2_1 = kron(permutedims(dGp),1.0I(size(D2,1)))*dD2;
          dT2_2 = kron(1.0I(size(dGp,2)), D2)*dGp_xx;
          dT2 = dT2_1 + dT2_2;
          return dT2
        end
        function jacT3(x)
          dT3_1 = kron(permutedims(Gp*dh_x), 1.0I(size(Ms,1)))*dMs;
          dT3_2 = kron(permutedims(dh_x), Ms)*dGp;
          dT3_3 = kron(1.0I(size(dh_x,2)), Ms*Gp)*dh_xx;
          dT3 = dT3_1 + dT3_2 + dT3_3;

          return dT3
        end
        dFc_xx = jacT1(x) + jacT2(x) + jacT3(x);

        # dFc_ux
        # # first term
        D1 = kron(permutedims(Gp*h), 1.0I(size(Ms,1)));
        dGph_u = Gp*dh_u;
        dGph_u_t = dGph_u;
        dD1_u = fdJ_kronAB(permutedims(Gp*h), 1.0I(size(Ms,1)), dGph_u_t, zeros(size(dMs,1),length(u)));
        out1 = kron(permutedims(dMs), 1.0I(size(D1,1)))*dD1_u;
        # # second term
        D2 = kron(permutedims(h), Ms)*dGp;
        dD2_u = fdJ_kronAB(permutedims(h), Ms, dh_u, zeros(size(dMs,1), length(u)));
        out2 = kron(permutedims(dGp), 1.0I(size(D2,1)))*dD2_u;
        # # third term
        out3 = kron(1.0I(size(dh_x,2)), Ms*Gp)*dh_ux;
        dFc_ux = out1 + out2 + out3;

        # dFc_xu
        out1 = kron(permutedims(Gp*dh_u), 1.0I(size(Ms,1)))*dMs;
        out2 = kron(permutedims(dh_u), Ms)*dGp;
        out3 = kron(1.0I(length(u)), Ms*Gp)*dh_xu;
        dFc_xu = out1 + out2 + out3;


        # dFc_uu
        dFc_uu = kron(1.0I(size(dh_u,2)), Ms*Gp)*dh_uu;

        return dFc_xx, dFc_ux, dFc_xu, dFc_uu
      end
      dFc_xx, dFc_ux, dFc_xu, dFc_uu = Fcdiff2(xv,uv);

      function hessqdd(x,u)
        T1 = kron(permutedims(Fc + Fu), 1.0I(size(M_inv,1)));
        dT1 = fdJ_kronAB(permutedims(Fc + Fu), 1.0I(size(M_inv,1)), (dFc + dFu), zeros(size(dM_inv)));
        T2 = dFc + dFu;
        dT2 = dFc_xx + dFu_xx;

        out1 = kron(permutedims(dM_inv), 1.0I(size(T1,1)))*dT1;
        out2 = kron(1.0I(size(dM_inv,2)), T1)*dMinv_xx;
        out3 = kron(permutedims(T2), 1.0I(size(M_inv,1)))*dM_inv;
        out4 = kron(1.0I(size(T2,2)), M_inv)*dT2;

        dqdd_xx = out1 + out2 + out3 + out4;
        return dqdd_xx
      end
      dqdd_xx = hessqdd(xv,uv);

      function cross_qdd_ux(x,u)
        T1 = kron(permutedims(Fc + Fu), 1.0I(size(M_inv,1)));
        dT1_u = fdJ_kronAB(permutedims(Fc + Fu), 1.0I(size(M_inv,1)), dFc_u + dFu_u, zeros(size(M_inv,1)^2, length(u)));

        T2 = dFc + dFu;
        dT2_u = dFc_ux + dFu_ux;

        out1 = kron(permutedims(dM_inv), 1.0I(size(T1,1)))*dT1_u;
        out2 = kron(1.0I(length(x)), M_inv)*dT2_u;
        dqdd_ux = out1+out2;
        return dqdd_ux
      end
      dqdd_ux = cross_qdd_ux(xv,uv);

      function cross_qdd_xu(x,u)
        out1 = kron(permutedims(dFc_u + dFu_u), 1.0I(size(M_inv,1)))*dM_inv;
        out2 = kron(1.0I(length(u)), M_inv)*(dFc_xu + dFu_xu);
        dqdd_xu = out1 + out2;
        return dqdd_xu
      end
      dqdd_xu = cross_qdd_xu(xv,uv);

      function hessqdd_u(x,u)
        dqdd_uu = kron(1.0I(length(u)), M_inv)*(dFc_uu + dFu_uu);
        return dqdd_uu
      end
      dqdd_uu = hessqdd_u(xv,uv);

      function hessDyn(xuVar)
        x = xuVar[1:14*nB]; u = reshape(xuVar[14*nB+1:end], (6,nB));

        y = qdd;
        dy = dqdd;

        d2y_xx = dqdd_xx;#  hessqdd(x,u);
        d2y_xu = dqdd_xu;# cross_qdd_xu(x,u);
        d2y_ux = dqdd_ux;# cross_qdd_ux(x,u);
        d2y_uu = dqdd_uu;# hessqdd_u(x,u);

        d2y = [d2y_xx d2y_ux; d2y_xu d2y_uu];

        out1 = kron(1.0I(length(xuVar)), permutedims(y))*d2y;
        out2 = permutedims(dy)*dy;
        out = 2*(out1 + out2); # hessian of dynamics cost = sum(y.*y)
        return out;
      end
      dqdd_xuVar = hessDyn(xu_var);

      function hess_nonlinear_constraints()
        zLen = length(z);

        hess_dyn = hessDyn(xu_var); # hessian of dynamics cost
        _, hess_qNorm = quatNormConstr_diff(xu_var); # hessian of quatNormConstr cost
        _, hess_jLoc = jointLocConstr_diff(xu_var); # hessian of jointLocConstr cost
        println()
        # println("typeof(hess_dyn) = ", typeof(hess_dyn))
        # println("typeof(hess_qNorm) = ", typeof(hess_qNorm))
        # println("typeof(hess_jLoc) = ", typeof(hess_jLoc))

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
  prob.x = inGuess; # Initial guess.
  status = solveProblem(prob);

  trim_x = prob.x[1:length(x0)];
  gam = prob.x[end];
  trim_u = prob.x[length(x0)+1:end-1];
  return trim_x, trim_u, gam
  # return trim_x, gam
end
## testing the trim function
# x0Orig, u0Orig = getXU_0(j);
# nB = length(j) + 1;
# ix= zeros(Integer,20*(nB-1)+1);
# # WITH U=0
# ix[14*(nB-1)+1:20*(nB-1)] .= 1;
# trim_U = u0Orig;
# iy = zeros(Integer,7*(nB-1));
#
# out = trimDiff2(j,g,ix=ix, iy=iy);
# trim_x, trim_u, gam = out;
# println("trim_x = ", trim_x)
# trim_X = [j[1].RB1.x;trim_x]
# # trim_U = [zeros(6) reshape(trim_u,(6,nB-1))];
# # println("trim_U = ", trim_U[:,2])
# # trim_U = zeros(6,nB);
# println("qconstr = ", norm(trim_x[4:7])-1)
# println("ẍ =", norm(fxdot(trim_X,trim_U, j, g)))
