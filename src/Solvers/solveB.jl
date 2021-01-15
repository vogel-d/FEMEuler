#Lösen der Akkustik-Gleichung mit Buoyancy
function solveB!(p::femProblem,
                 Fp::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}, Fv::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}, Fb::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},
                 tstart::Float64, dt::Float64, nt::Float64, method::Symbol)

  iterations=nt/dt;

  Spv=p.stiffM[:pv];
  Svp=p.stiffM[:vp];
  Svb=p.stiffM[:vb];
  Sbv=p.stiffM[:bv];

  yP=copy(p.solution[tstart].p);
  yV=copy(p.solution[tstart].v);
  yB=copy(p.solution[tstart].b);

  np=p.degFBoundary[p.femType[:p][1]].num
  nv=p.degFBoundary[p.femType[:v][1]].num
  nb=p.degFBoundary[p.femType[:b][1]].num

  if method==:euler
    for i in 1:iterations
      yV[1:nv]=yV[1:nv]+dt*(Fv\(-Svp*yP+Svb*yB));
      yP[1:np]=yP[1:np]+dt*(Fp\(-cs2*(Spv*yV)));
      yB[1:nb]=yB[1:nb]+dt*(Fb\(-N2*(Sbv*yV)));
    end

  elseif method==:rk4
    v1=copy(p.solution[tstart].v); v2=copy(p.solution[tstart].v); v3=copy(p.solution[tstart].v); v4=copy(p.solution[tstart].v);
    p1=copy(p.solution[tstart].p); p2=copy(p.solution[tstart].p); p3=copy(p.solution[tstart].p); p4=copy(p.solution[tstart].p);
    b1=copy(p.solution[tstart].b); b2=copy(p.solution[tstart].b); b3=copy(p.solution[tstart].b); b4=copy(p.solution[tstart].b);

    for i in 1:iterations
      v1[1:nv]=Fv\(-Svp*yP+Svb*yB);
      p1[1:np]=Fp\(-cs2*(Spv*yV));
      b1[1:nb]=Fb\(-N2*(Sbv*yV));

      v2[1:nv]=Fv\(-Svp*(yP+0.5*dt*p1)+Svb*(yB+0.5*dt*b1));
      p2[1:np]=Fp\(-cs2*(Spv*(yV+0.5*dt*v1)));
      b2[1:nb]=Fb\(-N2*(Sbv*(yV+0.5*dt*v1)));

      v3[1:nv]=Fv\(-Svp*(yP+0.5*dt*p2)+Svb*(yB+0.5*dt*b2));
      p3[1:np]=Fp\(-cs2*(Spv*(yV+0.5*dt*v2)));
      b3[1:nb]=Fb\(-N2*(Sbv*(yV+0.5*dt*v2)));

      v4[1:nv]=Fv\(-Svp*(yP+dt*p3)+Svb*(yB+dt*b3));
      p4[1:np]=Fp\(-cs2*(Spv*(yV+dt*v3)));
      b4[1:nb]=Fb\(-N2*(Sbv*(yV+dt*v3)));


      yV[1:nv]=yV[1:nv]+dt*(1/6)*(v1[1:nv]+2*v2[1:nv]+2*v3[1:nv]+v4[1:nv]);
      yP[1:np]=yP[1:np]+dt*(1/6)*(p1[1:np]+2*p2[1:np]+2*p3[1:np]+p4[1:np]);
      yB[1:nb]=yB[1:nb]+dt*(1/6)*(b1[1:nb]+2*b2[1:nb]+2*b3[1:nb]+b4[1:nb]);
    end
  else
    error("Keine zulässige Methode! Mögliche Methoden sind :euler und :rk4.");
  end

  p.solution[tstart+nt]=createSolution(yV,yP,yB);
end
