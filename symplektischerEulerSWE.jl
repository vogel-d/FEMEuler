function symplektischerEuler!(y::solution,p::femProblem,w::solution,ns::Int64,dtau::Float64,
                              velOld::Array{Float64,1},hVS::Array{Float64,1},hS::Array{Float64,1})
  gamma=0.01;
  numh=p.degFBoundary[p.femType[:h][1]].num
  numhV=p.degFBoundary[p.femType[:hV][1]].num

  for i in 1:ns
    copy!(velOld,y.hV);
    velocity!(hVS,p,y.h);
    @. y.hV[1:numhV]+=dtau*(hVS+w.hV[1:numhV]);
    height!(hS,p,(1+gamma)*y.hV-gamma*velOld);
    @. y.h[1:numh]+=dtau*(hS+w.h[1:numh]);
  end

  return nothing;
end

function velocity!(hVS::Array{Float64,1},p::femProblem, yh::Array{Float64,1})
  pres=projectPressureSWE(p.degFBoundary[p.femType[:p][1]],p.massMBoundary[p.femType[:p][1]],p.degFBoundary[p.femType[:h][1]],yh,p.mesh,p.kubPoints,p.kubWeights);
  hVS[:]=p.massM[p.femType[:hV][1]]\(-(p.stiffM[:vp]*pres));
  return nothing;
end

function height!(hS::Array{Float64,1},p::femProblem, yhV::Array{Float64,1})
  hS[:]=p.massM[p.femType[:h][1]]\(-p.stiffM[:h]*yhV);
  return nothing;
end
