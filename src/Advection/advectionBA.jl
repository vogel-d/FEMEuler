function advection(p::femProblem, gamma::Float64, Vfval::SparseVector{Float64,Int64}, Vfcomp::Symbol, cval::solution,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}})
  fTp=p.femType[:p]; fTb=p.femType[:b]; fTv=p.femType[:v];
  Fp=p.massM[fTp[1]]; Fb=p.massM[fTb[1]]; Fv=p.massM[fTv[1]];
  np=p.degFBoundary[fTp[1]].num; nv=p.degFBoundary[fTv[1]].num; nb=p.degFBoundary[fTb[1]].num;
  if p.taskRecovery
    if iszero(p.recoveryOrders[1])
      cR=recovery(p,1,fTp,cval.p);
    else
      cR=recovery(p,fTp,cval.p);
    end
    S=advectionStiff(p.degFBoundary[fTp[1]],nquadPhi[fTp[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTp[3]],nquadPhi[fTp[3]],cR,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.data);
    rCp=Fp\S;
    if iszero(p.recoveryOrders[2])
      cR=recovery(p,1,fTb,cval.b);
    else
      cR=recovery(p,fTb,cval.b);
    end
    S=advectionStiff(p.degFBoundary[fTb[1]],nquadPhi[fTb[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTb[3]],nquadPhi[fTb[3]],cR,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.data);
    rCb=Fb\S;
    if iszero(p.recoveryOrders[3])
      cR=recovery(p,2,fTv,cval.v);
    else
      cR=recovery(p,fTv,cval.v);
    end
    S=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.data);
    rCv=Fv\S;
  else
    S=advectionStiff(p.degFBoundary[fTp[1]],nquadPhi[fTp[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTp[1]],nquadPhi[fTp[1]],cval.p,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.data);
    rCp=Fp\S;
    S=advectionStiff(p.degFBoundary[fTb[1]],nquadPhi[fTb[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTb[1]],nquadPhi[fTb[1]],cval.b,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.data);

    rCb=Fb\S;
    S=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cval.v,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.data);

    rCv=Fv\S;
  end
  f=createSolution(length(cval.v),length(cval.p),length(cval.b));
  f.p[1:np]=rCp[:,1];
  f.b[1:nb]=rCb[:,1];
  f.v[1:nv]=rCv[:,1];
  return f;
end
