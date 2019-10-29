function splitExplicit(p::femProblem, gamma::Float64,Vfcomp::Symbol,
                      Vfval::SparseVector{Float64,Int64}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                      nquadPoints::Array{Array{Float64,2},1}, MIS::MIS,y0::solution,t0::Float64,dt::Float64,ns::Int64)
  stage=MIS.nStage;

  Y=Array{solution,1}(undef,stage+1);
  FY=Array{solution,1}(undef,stage);

  for i in 1:(stage+1)
    Y[i]=y0;
  end

  for i in 1:stage
    for j in 1:i
      Y[i+1]=MIS.alpha[i+1,j]*(Y[j]-y0)+Y[i+1];
    end

    FY[i]=advection(p,gamma,Vfval,Vfcomp,Y[i],nquadPoints,nquadPhi);

    fSlow=createSolution(length(y0.v),length(y0.p),length(y0.b));
    for j in 1:i
      fSlow=fSlow+MIS.beta[i+1,j]*FY[j]+(MIS.gamma[i+1,j]/dt)*(Y[j]-y0);
    end

    nsLoc=ceil(ns*MIS.d[i+1]);
    dtLoc=dt*MIS.d[i+1];
    dtauLoc=dtLoc/nsLoc;

    Y[i+1].p[p.degFBoundary[p.femType[:p][1]].num+1:end]=y0.p[p.degFBoundary[p.femType[:p][1]].num+1:end];
    Y[i+1].b[p.degFBoundary[p.femType[:b][1]].num+1:end]=y0.b[p.degFBoundary[p.femType[:b][1]].num+1:end];
    Y[i+1].v[p.degFBoundary[p.femType[:v][1]].num+1:end]=y0.v[p.degFBoundary[p.femType[:v][1]].num+1:end];

    symplektischerEuler!(Y[i+1],p,fSlow,nsLoc,dtauLoc);
  end

  return Y[stage+1];
end
