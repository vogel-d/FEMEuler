function symplektischerEuler!(y::solution,p::femProblem,w::solution,ns::Int64,dtau::Float64,
                              velOld::Array{Float64,1},rhoVS::Array{Float64,1},rhoS::Array{Float64,1},rhoThetaS::Array{Float64,1})
  gamma=0.01;
  numRho=p.degFBoundary[p.femType[:rho][1]].num
  numRhoV=p.degFBoundary[p.femType[:rhoV][1]].num
  numRhoTheta=p.degFBoundary[p.femType[:rhoTheta][1]].num

  for i in 1:ns
    copy!(velOld,y.rhoV);
    @. y.rhoV[1:numRhoV]+=dtau*(w.rhoV[1:numRhoV]);
    @. y.rho[1:numRho]+=dtau*(w.rho[1:numRho]);
    @. y.rhoTheta[1:numRhoTheta]+=dtau*(w.rhoTheta[1:numRhoTheta]);
  end

  return nothing;
end
