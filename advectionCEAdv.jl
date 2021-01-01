function advection(p::femProblem, gamma::Float64, y::solution, Vfval::SparseVector{Float64,Int64}, Vfcomp::Symbol,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                  MrT::SparseMatrixCSC{Float64,Int64}, MrV::SparseMatrixCSC{Float64,Int64})

  fTtheta=p.femType[:rhoTheta]; fTv=p.femType[:rhoV]; fTrho=p.femType[:rho];
  Fv=@views p.massM[fTv[1]]; Fth=@views p.massM[fTtheta[1]]; Frho=@views p.massM[fTrho[1]];
  nRhoV=p.degFBoundary[p.femType[:rhoV][1]].num;
  nRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num;
  nRho=p.degFBoundary[p.femType[:rho][1]].num;
  if p.advection
    cR=projectRhoChi(p,y.rho,y.rhoV,:rho,:rhoV,MrV);
    if p.taskRecovery
      if iszero(p.recoveryOrders[2])
        cR=recovery(p,2,fTv,cR);
        Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints,p.data.edgeData);

      else
        cR=recovery(p,fTv,cR);
        Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints,p.data.edgeData);
        #=
        cR=recovery(p,fTv,cR,2);
        Sv=advectionStiffR(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          fTv[3],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints,p.data.edgeData);
        =#
      end
      rCv=Fv\(Sv);
    else
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.data.edgeData);
      rCv=Fv\(Sv);
    end
    cR=projectRhoChi(p,y.rho,y.rhoTheta,:rho,:rhoTheta,MrT);
    if p.taskRecovery
      if iszero(p.recoveryOrders[1])
        cR=recovery(p,1,fTtheta,cR);
        Sth=advectionStiff(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                           p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                           p.degFBoundary[fTtheta[3]],nquadPhi[fTtheta[3]],cR,
                           gamma,p.mesh,p.kubPoints,p.kubWeights,
                           nquadPoints, p.data.edgeData);
      else
        cR=recovery(p,fTtheta,cR);
        Sth=advectionStiff(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                           p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                           p.degFBoundary[fTtheta[3]],nquadPhi[fTtheta[3]],cR,
                           gamma,p.mesh,p.kubPoints,p.kubWeights,
                           nquadPoints, p.data.edgeData);
        #=
        Sth=advectionStiffR(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                            p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                            fTtheta[3],cR,
                            gamma,p.mesh,p.kubPoints,p.kubWeights,
                            nquadPoints, p.data.edgeData);
        =#
      end
      rCth=Fth\(Sth);
    else
       Sth=advectionStiff(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints, p.data.edgeData);
      rCth=Fth\(Sth);
    end
    cR=y.rho
    if p.taskRecovery
      if iszero(p.recoveryOrders[3])
        cR=recovery(p,1,fTrho,cR);
        Srho=advectionStiff(p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],
                           p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                           p.degFBoundary[fTrho[3]],nquadPhi[fTrho[3]],cR,
                           gamma,p.mesh,p.kubPoints,p.kubWeights,
                           nquadPoints, p.data.edgeData);
      else
        cR=recovery(p,fTrho,cR);
        Srho=advectionStiff(p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],
                           p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                           p.degFBoundary[fTrho[3]],nquadPhi[fTrho[3]],cR,
                           gamma,p.mesh,p.kubPoints,p.kubWeights,
                           nquadPoints, p.data.edgeData);
        #=
        Srho=advectionStiffR(p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],
                             p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                             fTrho[3],cR,
                             gamma,p.mesh,p.kubPoints,p.kubWeights,
                             nquadPoints, p.data.edgeData);
        =#
      end
      rCrho=Frho\(Srho);
    else
       Srho=advectionStiff(p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints, p.data.edgeData);
      rCrho=Frho\(Srho);
    end
  else
    rCv=zeros(nRhoV);
    rCth=zeros(nRhoTheta);
    rCrho=zeros(nRho);
  end
  f=createSolution(length(y.rho),length(y.rhoV),length(y.rhoTheta),length(y.v),length(y.theta));
  @views f.rhoV[1:nRhoV]=rCv[:,1];
  @views f.rhoTheta[1:nRhoTheta]=rCth[:,1];
  @views f.rho[1:nRho]=rCrho[:,1];
  return f;
end
