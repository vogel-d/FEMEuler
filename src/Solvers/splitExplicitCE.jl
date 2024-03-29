function splitExplicit(y0::solution,Y::Array{solution,1},FY::Array{solution,1},SthY::Array{SparseMatrixCSC{Float64,Int64},1},
                      p::femProblem, gamma, #Sth::SparseMatrixCSC{Float64,Int64},
                      nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                      nquadPoints::Array{Array{Float64,2},1},MrT::SparseMatrixCSC{Float64,Int64}, MrV::SparseMatrixCSC{Float64,Int64},
                      MIS::MIS,t0::Float64,dt::Float64,ns::Int64)
  stage=MIS.nStage;
  #Y=Array{solution,1}(undef,stage+1);
  #FY=Array{solution,1}(undef,stage);
  #SthY=Array{SparseMatrixCSC{Float64,Int64},1}(undef,stage);

  numRho=p.degFBoundary[p.femType[:rho][1]].num
  numRhoV=p.degFBoundary[p.femType[:rhoV][1]].num
  numRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num

  velOld=Array{Float64,1}(undef,numRhoV)
  rhoVS=Array{Float64,1}(undef,numRhoV)
  rhoS=Array{Float64,1}(undef,numRho)
  rhoThetaS=Array{Float64,1}(undef,numRhoTheta)

  for i in 1:(stage+1)
    Y[i]=y0;
  end

  for i in 1:stage
    for j in 1:i
      Y[i+1]=MIS.alpha[i+1,j]*(Y[j]-y0)+Y[i+1];
    end
    FY[i], SthY[i]=advection(p,gamma,Y[i],nquadPoints,nquadPhi,MrT,MrV);
    #FY[i]=advection(p,gamma,Y[i],nquadPoints,nquadPhi,MrT,MrV);
    #SthY[i]=Sth;

    fSlow=createSolution(length(y0.rho),length(y0.rhoV),length(y0.rhoTheta),length(y0.v),length(y0.theta));
    SthSlow=spzeros(size(SthY[i],1),size(SthY[i],2))

    for j in 1:i
      fSlow=fSlow+MIS.beta[i+1,j]*FY[j]+(MIS.gamma[i+1,j]/dt)*(Y[j]-y0);
      SthSlow=SthSlow+MIS.beta[i+1,j]*SthY[j];
    end

    nsLoc=ceil(Int64,ns*MIS.d[i+1]);
    dtLoc=dt*MIS.d[i+1];
    dtauLoc=dtLoc/nsLoc;

    Y[i+1].rhoV[numRhoV+1:end]=y0.rhoV[numRhoV+1:end];
    Y[i+1].rho[numRho+1:end]=y0.rho[numRho+1:end];
    Y[i+1].rhoTheta[numRhoTheta+1:end]=y0.rhoTheta[numRhoTheta+1:end];

    symplektischerEuler!(Y[i+1],p,fSlow,SthSlow,nsLoc,dtauLoc,velOld,rhoVS,rhoS,rhoThetaS);
  end
  return Y[stage+1];
end
