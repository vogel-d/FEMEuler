function advection(p::femProblem, gamma, y::solution,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                  MrV::SparseMatrixCSC{Float64,Int64})

  fTv=p.femType[:hV];
  Fv=@views p.massM[fTv[1]];
  nhV=p.degFBoundary[p.femType[:hV][1]].num;

  if p.advection
    cR=projectRhoChi(p,y.h,y.hV,:h,:hV,MrV);
    if p.taskRecovery
      if iszero(p.recoveryOrders[2])
        cR=recovery(p,2,fTv,cR);
        Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                          p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.hV),
                          p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints,p.data);

      else
        cR=recovery(p,fTv,cR);
        Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                          p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.hV),
                          p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints,p.data);
        #=
        cR=recovery(p,fTv,cR,2);
        Sv=advectionStiffR(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                          p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.hV),
                          fTv[3],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints,p.data.edgeData);
        =#
      end
    else
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.hV),
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.data);
    end
    rCv=Fv\(Sv+p.stiffM[:fv]*y.hV);
  else
    rCv=zeros(nRhoV);
  end

  f=createSolution(length(y.h),length(y.hV),length(y.v));
  @views f.hV[1:nhV]=rCv[:,1];

  return f;
end
