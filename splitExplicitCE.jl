function splitExplicit(y0::solution,Y::Array{solution,1},FY::Array{solution,1},SthY::Array{SparseMatrixCSC{AbstractFloat,Int},1},
                      p::femProblem, gamma::AbstractFloat,
                      nquadPhi::Dict{Symbol, Array{Array{Array{AbstractFloat,1},2},1}},
                      nquadPoints::Array{Array{AbstractFloat,2},1},MrT::SparseMatrixCSC{AbstractFloat,Int}, MrV::SparseMatrixCSC{AbstractFloat,Int},
                      MIS::MIS,t0::AbstractFloat,dt::AbstractFloat,ns::Int)
  stage=MIS.nStage;
  #Y=Array{solution,1}(undef,stage+1);
  #FY=Array{solution,1}(undef,stage);
  #SthY=Array{SparseMatrixCSC{AbstractFloat,Int},1}(undef,stage);

  numRho=p.degFBoundary[p.femType[:rho][1]].num
  numRhoV=p.degFBoundary[p.femType[:rhoV][1]].num
  numRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num

  velOld=Array{AbstractFloat,1}(undef,numRhoV)
  rhoVS=Array{AbstractFloat,1}(undef,numRhoV)
  rhoS=Array{AbstractFloat,1}(undef,numRho)
  rhoThetaS=Array{AbstractFloat,1}(undef,numRhoTheta)

  for i in 1:(stage+1)
    Y[i]=y0;
  end

  for i in 1:stage
    for j in 1:i
      Y[i+1]=MIS.alpha[i+1,j]*(Y[j]-y0)+Y[i+1];
    end

    FY[i], SthY[i]=advection(p,gamma,Y[i],nquadPoints,nquadPhi,MrT,MrV);

    fSlow=createSolution(length(y0.rho),length(y0.rhoV),length(y0.rhoTheta),length(y0.v),length(y0.theta));
    SthSlow=spzeros(size(SthY[i],1),size(SthY[i],2))

    for j in 1:i
      fSlow=fSlow+MIS.beta[i+1,j]*FY[j]+(MIS.gamma[i+1,j]/dt)*(Y[j]-y0);
      SthSlow=SthSlow+MIS.beta[i+1,j]*SthY[j];
    end

    nsLoc=ceil(Int,ns*MIS.d[i+1]);
    dtLoc=dt*MIS.d[i+1];
    dtauLoc=dtLoc/nsLoc;

    Y[i+1].rhoV[numRhoV+1:end]=y0.rhoV[numRhoV+1:end];
    Y[i+1].rho[numRho+1:end]=y0.rho[numRho+1:end];
    Y[i+1].rhoTheta[numRhoTheta+1:end]=y0.rhoTheta[numRhoTheta+1:end];

    symplektischerEuler!(Y[i+1],p,fSlow,SthSlow,nsLoc,dtauLoc,velOld,rhoVS,rhoS,rhoThetaS);
  end
  return Y[stage+1];
end
