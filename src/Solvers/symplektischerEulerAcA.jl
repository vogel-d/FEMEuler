function symplektischerEuler!(y::solution,p::femProblem, w::solution,ns::Float64,dtau::Float64,
                              velOld::Array{Float64,1},vS::Array{Float64,1},pS::Array{Float64,1},bS::Array{Float64,1})

  gamma=0.05;

  np=p.degFBoundary[p.femType[:p][1]].num
  nv=p.degFBoundary[p.femType[:v][1]].num

  for i in 1:ns
    copy!(velOld,y.v);
    velocity!(vS,p,y.p)
    @. y.v[1:nv]+=dtau*(vS+w.v[1:nv]);
    pressure!(pS,p,(1+gamma)*y.v-gamma*velOld);
    @. y.p[1:np]+=dtau*(pS+w.p[1:np]);

  end

  return nothing;
end

function velocity!(vS::Array{Float64,1},p::femProblem, yP::Array{Float64,1})
  Fv=p.massM[p.femType[:v][1]];
  Svp=p.stiffM[:vp];
  vS[:]=Fv\(-Svp*yP)
  return nothing;
end

function pressure!(pS::Array{Float64,1}, p::femProblem, yV::Array{Float64,1})
  Fp=p.massM[p.femType[:p][1]];
  Spv=p.stiffM[:pv];
  pS[:]=Fp\(-cs2*(Spv*yV));
  return nothing;
end
