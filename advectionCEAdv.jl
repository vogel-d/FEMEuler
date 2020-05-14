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
      cR=recovery(p,fTv,cR,p.stencil);
      #vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][3]],cR,p.femType[:rhoV][3],"testRecoveryV")
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                        p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.edgeData);
      rCv=Fv\(Sv);
    else
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.edgeData);
      rCv=Fv\(Sv);
    end
    cR=projectRhoChi(p,y.rho,y.rhoTheta,:rho,:rhoTheta,MrT);
    if p.taskRecovery
      #vtk(p.mesh,p.degFBoundary[p.femType[:rhoTheta][1]],cR,p.femType[:rhoTheta][1],"testSphereTh")
      cR=recovery(p,fTtheta,cR);
      #vtkRecovery(p.mesh,2,p.degFBoundary[p.femType[:rhoTheta][3]],cR,p.femType[:rhoTheta][3],"testSphereRecoveryTh")
      Sth=advectionStiffR(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                         fTtheta[3],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.edgeData);
      rCth=Fth\(Sth);
    else
       Sth=advectionStiff(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints, p.edgeData);
      rCth=Fth\(Sth);
      #vtk(p.mesh,p.degFBoundary[p.femType[:rhoTheta][1]],Sth,p.femType[:rhoTheta][1],"testSphereThAdvS$zvtk")
      #vtk(p.mesh,p.degFBoundary[p.femType[:rhoTheta][1]],rCth,p.femType[:rhoTheta][1],"testSphereThAdvrCth$zvtk")
    end
    cR=y.rho
    if p.taskRecovery
      cR=recovery(p,fTrho,cR);
      #vtkRecovery(p.mesh,2,p.degFBoundary[p.femType[:rho][3]],cR,p.femType[:rho][3],"testRecoveryRho")
      Srho=advectionStiffR(p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],
                         p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                         fTrho[3],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.edgeData);
      rCrho=Frho\(Srho);
    else
       Srho=advectionStiff(p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],
                          p.degFBoundary[Vfcomp],nquadPhi[Vfcomp],Vfval,
                          p.degFBoundary[fTrho[1]],nquadPhi[fTrho[1]],cR,
                          gamma,p.mesh,p.kubPoints,p.kubWeights,
                          nquadPoints, p.edgeData);
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
