#Lösen der Akkustik-Gleichung mit Buoyancy
function RKadvection!(p::femProblem, gamma::Float64, Vfcomp::Symbol,
                      Vfval::SparseVector{Float64,Int64}, nquadPhi::Dict{Symbol, Array{Array{Array{Float64,1},2},1}},
                      nquadPoints::Array{Array{Float64,2},1}, Time::Float64, dt::Float64, EndTime::Float64,
                      method::Symbol)

  timesteps=(EndTime-Time)/dt;
  !isequal(mod(timesteps,1),0) && @warn("Es existiert kein ganzzahliges k, sodass k*dt=$EndTime. Endzeitpunkt dieser Berechnung: $(Time + floor(timesteps)*dt)")

  timesteps=Int(floor(timesteps));

  Fp=p.massM[p.femType[:p][1]];
  Fv=p.massM[p.femType[:v][1]];
  Fb=p.massM[p.femType[:b][1]];

  #f=copy(p.solution[Time])

  np=p.degFBoundary[p.femType[:p][1]].num
  nv=p.degFBoundary[p.femType[:v][1]].num
  nb=p.degFBoundary[p.femType[:b][1]].num

  if method==:euler
    f=createSolution(length(p.solution[Time].v),length(p.solution[Time].p),length(p.solution[Time].b));
    f.p=copy(p.solution[time].p);
    f.b=copy(p.solution[time].b);
    f.v=copy(p.solution[time].v);

    for i in 1:timesteps
      Sp, Sb, Sv = advection(p,gamma,Vfval,Vfcomp,f,nquadPoints,nquadPhi);
      f.v[1:nv]=f.v[1:nv]+dt*(Fv\Sp);
      f.p[1:np]=f.p[1:np]+dt*(Fp\Sb);
      f.b[1:nb]=f.b[1:nb]+dt*(Fb\Sv);

      p.solution[Time+dt]=copy(f);
      Time+=dt;
    end
  elseif method==:rk4
    f =createSolution(length(p.solution[Time].v),length(p.solution[Time].p),length(p.solution[Time].b));
    f1=createSolution(length(p.solution[Time].v),length(p.solution[Time].p),length(p.solution[Time].b));
    f2=createSolution(length(p.solution[Time].v),length(p.solution[Time].p),length(p.solution[Time].b));
    f3=createSolution(length(p.solution[Time].v),length(p.solution[Time].p),length(p.solution[Time].b));
    f4=createSolution(length(p.solution[Time].v),length(p.solution[Time].p),length(p.solution[Time].b));

    f.p=copy(p.solution[Time].p);
    f.b=copy(p.solution[Time].b);
    f.v=copy(p.solution[Time].v);

    for i in 1:timesteps
      Sp, Sb, Sv = advection(p,gamma,Vfval,Vfcomp,f,nquadPoints,nquadPhi);
      f1.v[1:nv]=Fv\Sv
      f1.p[1:np]=Fp\Sp
      f1.b[1:nb]=Fb\Sb

      Sp, Sb, Sv = advection(p,gamma,Vfval,Vfcomp,f+0.5*dt*f1,nquadPoints,nquadPhi);
      f2.v[1:nv]=Fv\Sv
      f2.p[1:np]=Fp\Sp
      f2.b[1:nb]=Fb\Sb

      Sp, Sb, Sv = advection(p,gamma,Vfval,Vfcomp,f+0.5*dt*f2,nquadPoints,nquadPhi);
      f3.v[1:nv]=Fv\Sv
      f3.p[1:np]=Fp\Sp
      f3.b[1:nb]=Fb\Sb

      Sp, Sb, Sv = advection(p,gamma,Vfval,Vfcomp,f+dt*f3,nquadPoints,nquadPhi);
      f4.v[1:nv]=Fv\Sv
      f4.p[1:np]=Fp\Sp
      f4.b[1:nb]=Fb\Sb

      f.v[1:nv]=f.v[1:nv]+dt*(1/6)*(f1.v[1:nv]+2*f2.v[1:nv]+2*f3.v[1:nv]+f4.v[1:nv]);
      f.p[1:np]=f.p[1:np]+dt*(1/6)*(f1.p[1:np]+2*f2.p[1:np]+2*f3.p[1:np]+f4.p[1:np]);
      f.b[1:nb]=f.b[1:nb]+dt*(1/6)*(f1.b[1:nb]+2*f2.b[1:nb]+2*f3.b[1:nb]+f4.b[1:nb]);



      p.solution[Time+dt]=copy(f);
      Time+=dt;
    end
  else
    error("Keine zulässige Methode! Mögliche Methoden sind :euler und :rk4.");
  end
end
