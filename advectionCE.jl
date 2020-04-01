function advection(p::femProblem, gamma::Float64, y::solution,
                  nquadPoints::Array{Array{Float64,2},1}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                  MrT::SparseMatrixCSC{Float64,Int64}, MrV::SparseMatrixCSC{Float64,Int64})

  fTtheta=p.femType[:rhoTheta]; fTv=p.femType[:rhoV];
  Fv=@views p.massM[fTv[1]];
  nRhoV=p.degFBoundary[p.femType[:rhoV][1]].num;
  nRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num;

  if p.advection
    cR=projectRhoChi(p,y.rho,y.rhoV,:rho,:rhoV,MrV);
    if p.taskRecovery
      cR=recovery(p,2,fTv,cR);
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.rhoV),
                        p.degFBoundary[fTv[3]],nquadPhi[fTv[3]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.edgeData);
      rCv=Fv\(Sv+p.stiffM[:fv]*y.rhoV);
      #=
      pS=false;
      vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],Sv,p.femType[:rhoV][1],"advA", printSpherical=pS)
      resa=projectRhoChi(p,p.solution[0.0].rho,Sv,:rho,:rhoV,MrV)
      vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],resa,p.femType[:rhoV][1],"advARho", printSpherical=pS)
      vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],p.stiffM[:fv]*y.rhoV,p.femType[:rhoV][1],"coriolisA", printSpherical=pS)
      resc=projectRhoChi(p,p.solution[0.0].rho,p.stiffM[:fv]*p.solution[0.0].rhoV,:rho,:rhoV,MrV)
      vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],resc,p.femType[:rhoV][1],"coriolisRho", printSpherical=pS)
      #vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],Sv-p.stiffM[:fv]*y.rhoV,p.femType[:rhoV][1],"advCoriolisA", printSpherical=true)
      #vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],resa-resc,"advCoriolisA", printSpherical=true)
      resca=projectRhoChi(p,p.solution[0.0].rho,-Sv-p.stiffM[:fv]*y.rhoV,:rho,:rhoV,MrV)
      vtk(p.mesh,p.degFBoundary[p.femType[:rhoV][1]],resca,p.femType[:rhoV][1],"advCoriolisARho", printSpherical=pS)
      println(gfcjzd)
      =#
    else
      Sv=advectionStiff(p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],sparse(y.rhoV),
                        p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],cR,
                        gamma,p.mesh,p.kubPoints,p.kubWeights,
                        nquadPoints,p.edgeData);
      rCv=Fv\(Sv+p.stiffM[:fv]*y.rhoV);
    end
    cR=projectRhoChi(p,y.rho,y.rhoTheta,:rho,:rhoTheta,MrT);
    if p.taskRecovery
      cR=recovery(p,1,fTtheta,cR,:theta);
      Sth=advectionStiffMatrix(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],y.rhoV,
                         p.degFBoundary[fTtheta[3]],nquadPhi[fTtheta[3]],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.edgeData);
    else
      Sth=advectionStiffMatrix(p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],
                         p.degFBoundary[fTv[1]],nquadPhi[fTv[1]],y.rhoV,
                         p.degFBoundary[fTtheta[1]],nquadPhi[fTtheta[1]],cR,
                         gamma,p.mesh,p.kubPoints,p.kubWeights,
                         nquadPoints, p.edgeData);
    end
  else
    rCv=zeros(nRhoV);
    Sth=spzeros(nRhoTheta,length(y.rhoV));
  end

  f=createSolution(length(y.rho),length(y.rhoV),length(y.rhoTheta),length(y.v),length(y.theta));
  @views f.rhoV[1:nRhoV]=rCv[:,1];

  return f, Sth;
end
