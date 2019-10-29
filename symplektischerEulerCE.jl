function symplektischerEuler!(y::solution,p::femProblem,w::solution,Sth::SparseMatrixCSC{Float64,Int64},ns::Int64,dtau::Float64)
  gamma=0.01;
  numRho=p.degFBoundary[p.femType[:rho][1]].num
  numRhoV=p.degFBoundary[p.femType[:rhoV][1]].num
  numRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num

  for i in 1:ns
    velOld=copy(y.rhoV);
    rhoVS=velocity(p,y.rho,y.rhoTheta);
    y.rhoV[1:numRhoV]=y.rhoV[1:numRhoV]+dtau*(rhoVS+w.rhoV[1:numRhoV]);
    rhoS,rhoThetaS=rhoTheta(p,(1+gamma)*y.rhoV-gamma*velOld, Sth);
    y.rho[1:numRho]=y.rho[1:numRho]+dtau*(rhoS+w.rho[1:numRho]);
    y.rhoTheta[1:numRhoTheta]=y.rhoTheta[1:numRhoTheta]+dtau*(rhoThetaS+w.rhoTheta[1:numRhoTheta]);
  end

  return nothing;
end

function velocity(p::femProblem, yRho::Array{Float64,1},yRhoTheta::Array{Float64,1})
  pres=projectPressure(p,p.femType[:p][1],p.femType[:rhoTheta][1],yRhoTheta);
  return p.massM[p.femType[:rhoV][1]]\(-(p.stiffM[:vp]*pres+9.81*(p.stiffM[:vrho]*yRho)));
end

function rhoTheta(p::femProblem, yRhoV::Array{Float64,1}, STh::SparseMatrixCSC{Float64,Int64})
  rhoS=p.massM[p.femType[:rho][1]]\(-p.stiffM[:rho]*yRhoV);
  rhoThetaS=p.massM[p.femType[:rhoTheta][1]]\(STh*yRhoV);
  return rhoS, rhoThetaS;
end
