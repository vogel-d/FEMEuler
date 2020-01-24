function symplektischerEuler!(y::solution,p::femProblem, w::solution,ns::AbstractFloat,dtau::AbstractFloat,
                              velOld::Array{AbstractFloat,1},vS::Array{AbstractFloat,1},pS::Array{AbstractFloat,1},bS::Array{AbstractFloat,1})

  gamma=0.05;

  np=p.degFBoundary[p.femType[:p][1]].num
  nv=p.degFBoundary[p.femType[:v][1]].num
  nb=p.degFBoundary[p.femType[:b][1]].num

  for i in 1:ns
    copy!(velOld,y.v);
    velocity!(vS,p,y.p,y.b)
    @. y.v[1:nv]+=dtau*(vS+w.v[1:nv]);
    pressureBuoyancy!(pS,bS,p,(1+gamma)*y.v-gamma*velOld);
    @. y.p[1:np]+=dtau*(pS+w.p[1:np]);
    @. y.b[1:nb]+=dtau*(bS+w.b[1:nb]);

  end

  return nothing;
end

function velocity!(vS::Array{AbstractFloat,1},p::femProblem, yP::Array{AbstractFloat,1},yB::Array{AbstractFloat,1})
  Fv=p.massM[p.femType[:v][1]];
  Svp=p.stiffM[:vp];
  Svb=p.stiffM[:vb];
  vS[:]=Fv\(-Svp*yP+Svb*yB)
  return nothing;
end

function pressureBuoyancy!(pS::Array{AbstractFloat,1},bS::Array{AbstractFloat,1},p::femProblem, yV::Array{AbstractFloat,1})
  cs2=115600;
  N2=1.e-4;

  Fp=p.massM[p.femType[:p][1]];
  Fb=p.massM[p.femType[:b][1]];

  Spv=p.stiffM[:pv];
  Sbv=p.stiffM[:bv];

  pS[:]=Fp\(-cs2*(Spv*yV));
  bS[:]=Fb\(-N2*(Sbv*yV));

  return nothing;
end
