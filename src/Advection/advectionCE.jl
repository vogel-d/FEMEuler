function advection(p::femProblem, gamma, y::solution,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                  MrT::SparseMatrixCSC{Float64,Int64}, MrV::SparseMatrixCSC{Float64,Int64})

  fTtheta=p.femType[:rhoTheta]; fTv=p.femType[:rhoV];
  Fv=@views p.massM[fTv[1]];
  nRhoV=p.degFBoundary[p.femType[:rhoV][1]].num;
  nRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num;
  if p.advection
    cR=projectRhoChi(p,y.rho,y.rhoV,:rho,:rhoV,MrV);
    if p.taskRecovery
      if iszero(p.recoveryOrders[2])
        cR=recovery(p,2,fTv,cR);
      else
        cR=recovery(p,fTv,cR);
      end
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.rhoV),
                        p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.data);
    else
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.rhoV),
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.data);
    end
    if p.type==:shallow
      rCv=Fv\(Sv+p.stiffM[:fv]*y.rhoV);
    else
      rCv=Fv\(Sv);
    end
    cR=projectRhoChi(p,y.rho,y.rhoTheta,:rho,:rhoTheta,MrT);
    if p.taskRecovery
      if iszero(p.recoveryOrders[1])
        cR=recovery(p,1,fTtheta,cR,:theta);
      else
        cR=recovery(p,fTtheta,cR);
      end
      Sth=advectionStiffMatrix(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],y.rhoV,
                         p.degFBoundary[fTtheta[3]],nquadPhi[fTtheta[3]],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.data);
    else
      Sth=advectionStiffMatrix(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],y.rhoV,
                         p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.data);
    end
  else
    rCv=zeros(nRhoV);
    Sth=spzeros(nRhoTheta,length(y.rhoV));
  end
  f=createSolution(length(y.rho),length(y.rhoV),length(y.rhoTheta),length(y.v),length(y.theta));
  @views f.rhoV[1:nRhoV]=rCv[:,1];

  return f, Sth;
end
