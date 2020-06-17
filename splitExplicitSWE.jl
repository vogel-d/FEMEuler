function splitExplicit(y0::solution,Y::Array{solution,1},FY::Array{solution,1},
                      p::femProblem, gamma,
                      nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                      nquadPoints::Array{Array{Float64,2},1}, MrV::SparseMatrixCSC{Float64,Int64},
                      MIS::MIS,t0::Float64,dt::Float64,ns::Int64)
  stage=MIS.nStage;
  #Y=Array{solution,1}(undef,stage+1);
  #FY=Array{solution,1}(undef,stage);

  numh=p.degFBoundary[p.femType[:h][1]].num
  numhV=p.degFBoundary[p.femType[:hV][1]].num

  velOld=Array{Float64,1}(undef,numhV)
  hVS=Array{Float64,1}(undef,numhV)
  hS=Array{Float64,1}(undef,numh)

  for i in 1:(stage+1)
    Y[i]=y0;
  end

  for i in 1:stage
    for j in 1:i
      Y[i+1]=MIS.alpha[i+1,j]*(Y[j]-y0)+Y[i+1];
    end

    FY[i]=advection(p,gamma,Y[i],nquadPoints,nquadPhi,MrV);

    fSlow=createSolution(length(y0.h),length(y0.hV), length(y0.v));

    for j in 1:i
      fSlow=fSlow+MIS.beta[i+1,j]*FY[j]+(MIS.gamma[i+1,j]/dt)*(Y[j]-y0);
    end

    nsLoc=ceil(Int64,ns*MIS.d[i+1]);
    dtLoc=dt*MIS.d[i+1];
    dtauLoc=dtLoc/nsLoc;

    Y[i+1].hV[numhV+1:end]=y0.hV[numhV+1:end];
    Y[i+1].h[numh+1:end]=y0.h[numh+1:end];

    symplektischerEuler!(Y[i+1],p,fSlow,nsLoc,dtauLoc,velOld,hVS,hS);
  end
  return Y[stage+1];
end
