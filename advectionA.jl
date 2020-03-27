function advection(p::femProblem, gamma::Float64, Vfval::SparseVector{Float64,Int64}, Vfcomp::Symbol, cval::solution,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}})
  fTp=p.femType[:p]; fTb=p.femType[:b]; fTv=p.femType[:v];
  Fp=p.massM[fTp[1]]; Fb=p.massM[fTb[1]]; Fv=p.massM[fTv[1]];
  np=p.degFBoundary[fTp[1]].num; nv=p.degFBoundary[fTv[1]].num; nb=p.degFBoundary[fTb[1]].num;
  if p.taskRecovery
    cR=recovery(p,1,fTp,cval.p);
    Sp=advectionStiff(p.degFBoundary[fTp[1]],nquadPhi[fTp[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTp[3]],nquadPhi[fTp[3]],cR,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.edgeData);

    cR=recovery(p,1,fTb,cval.b);
    Sb=advectionStiff(p.degFBoundary[fTb[1]],nquadPhi[fTb[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTb[3]],nquadPhi[fTb[3]],cR,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.edgeData);

    cR=recovery(p,2,fTv,cval.v);
    Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.edgeData);
  else
    Sp=advectionStiff(p.degFBoundary[fTp[1]],nquadPhi[fTp[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTp[1]],nquadPhi[fTp[1]],cval.p,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.edgeData);

    Sb=advectionStiff(p.degFBoundary[fTb[1]],nquadPhi[fTb[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTb[1]],nquadPhi[fTb[1]],cval.b,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.edgeData);

    Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                     p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                     p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cval.v,
                     gamma,p.mesh,p.kubPoints,p.kubWeights,
                     nquadPoints,p.edgeData);
  end
  return Sp,Sb,Sv
end
