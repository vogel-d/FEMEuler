function symplektischerEuler!(y::solution,p::femProblem,w::solution,Sth::SparseMatrixCSC{Float64,Int64},ns::Int64,dtau::Float64,
                              velOld::Array{Float64,1},rhoVS::Array{Float64,1},rhoS::Array{Float64,1},rhoThetaS::Array{Float64,1})
  gamma=0.01;
  numRho=p.degFBoundary[p.femType[:rho][1]].num
  numRhoV=p.degFBoundary[p.femType[:rhoV][1]].num
  numRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num

  for i in 1:ns
    copy!(velOld,y.rhoV);
    velocity!(rhoVS,p,y.rho,y.rhoTheta);
    @. y.rhoV[1:numRhoV]+=dtau*(rhoVS+w.rhoV[1:numRhoV]);
    rhoTheta!(rhoS,rhoThetaS,p,(1+gamma)*y.rhoV-gamma*velOld, Sth);
    @. y.rho[1:numRho]+=dtau*(rhoS+w.rho[1:numRho]);
    @. y.rhoTheta[1:numRhoTheta]+=dtau*(rhoThetaS+w.rhoTheta[1:numRhoTheta]);
  end

  return nothing;
end

function velocity!(rhoVS::Array{Float64,1},p::femProblem, yRho::Array{Float64,1},yRhoTheta::Array{Float64,1})
  pres=projectPressure(p.degFBoundary[p.femType[:p][1]],p.massMBoundary[p.femType[:p][1]],p.degFBoundary[p.femType[:rhoTheta][1]],yRhoTheta,p.mesh,p.kubPoints,p.kubWeights,Val(p.type));
  rhoVS[:]=p.massM[p.femType[:rhoV][1]]\(-(p.stiffM[:vp]*pres+9.81*(p.stiffM[:vrho]*yRho)));
  return nothing;
end

function rhoTheta!(rhoS::Array{Float64,1},rhoThetaS::Array{Float64,1},p::femProblem, yRhoV::Array{Float64,1}, STh::SparseMatrixCSC{Float64,Int64})
  rhoS[:]=p.massM[p.femType[:rho][1]]\(-p.stiffM[:rho]*yRhoV);
  rhoThetaS[:]=p.massM[p.femType[:rhoTheta][1]]\(STh*yRhoV);
  return nothing;
end
