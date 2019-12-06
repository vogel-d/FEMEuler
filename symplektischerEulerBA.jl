function symplektischerEuler!(y::solution,p::femProblem, w::solution,ns::Float64,dtau::Float64)

  gamma=0.05;

  np=p.degFBoundary[p.femType[:p][1]].num
  nv=p.degFBoundary[p.femType[:v][1]].num
  nb=p.degFBoundary[p.femType[:b][1]].num

  for i in 1:ns
    velOld=copy(y.v);
    y.v[1:nv]=y.v[1:nv]+dtau*(w.v[1:nv]);
    y.p[1:np]=y.p[1:np]+dtau*(w.p[1:np]);
    y.b[1:nb]=y.b[1:nb]+dtau*(w.b[1:nb]);
    #=
    velOld=copy(y.v);
    y.v[1:nv]=y.v[1:nv]+dtau*(velocity(p,y.p,y.b)+w.v[1:nv]);
    pS,bS=pressureBuoyancy(p,(1+gamma)*y.v-gamma*velOld);
    y.p[1:np]=y.p[1:np]+dtau*(pS+w.p[1:np]);
    y.b[1:nb]=y.b[1:nb]+dtau*(bS+w.b[1:nb]);
    =#
  end

  return nothing;
end

function velocity(p::femProblem, yP::Array{Float64,1},yB::Array{Float64,1})
  Fv=p.massM[p.femType[:v][1]];
  Svp=p.stiffM[:vp];
  Svb=p.stiffM[:vb];

  return Fv\(-Svp*yP+Svb*yB);
end

function pressureBuoyancy(p::femProblem, yV::Array{Float64,1})
  cs2=115600;
  N2=1.e-4;

  Fp=p.massM[p.femType[:p][1]];
  Fb=p.massM[p.femType[:b][1]];

  Spv=p.stiffM[:pv];
  Sbv=p.stiffM[:bv];

  yP=Fp\(-cs2*(Spv*yV));
  yB=Fb\(-N2*(Sbv*yV));

  return yP, yB;
end
