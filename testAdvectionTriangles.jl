include("modulesBA.jl")
include("symplektischerEulerA.jl")

function testAdvectionTri()
  filename = "testosw"

  #order: comp, compHigh, compRec, compDG
  femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG0, :P1, :DG1, :DG0]);
  #femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:P1, :P1, :DG1, :DG1]);
  #femType=Dict(:p=>[:DG1, :P1, :DG1, :DG0], :v=>[:RT1, :VecP1, :VecDG1, :RT0B], :b=>[:DG1, :P1, :DG1, :DG0]);
  Vfcomp=:RT0
  #Vfcomp=:RT1

  taskRecovery=false;

  #m=generateTriMesh(300,20,:periodic,:periodic,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
  #m=generateTriMesh(100,100,:periodic,:periodic,0.0,20000.0,0.0,20000.0); #(east/west, top/bottom)
  m=generateTriMesh(10,10,:periodic,:periodic); #(east/west, top/bottom)
  #adaptGeometry!(m,(0.3,0.3),false); #sin perbutation
  p=femProblem(m, femType, t=:boussinesq, taskRecovery=taskRecovery);

  gamma=0.0;
  UMax=1.0;
  #UMax=20.0; #UMax determines the advection in x direction
  #UMax=0.0;
  MISMethod=MIS(:MIS_Euler);

  #dt=1000.0;
  dt=0.5;
  ns=19;
  #EndTime=20*dt;
  EndTime=300.0;
  nIter=Int64(EndTime/dt);

  #start function
  xR=m.geometry.r[1]; xL=m.geometry.l[1]; yR=m.geometry.r[2]; yL=m.geometry.l[2]
  b0=0.01; H=10000; A=5000;
  xM=0.5*(xL+xR);
  fb1(x,y)=b0*sin(pi*y/H)/(1+((x-xM)/A)^2);
  fb2(x,y)=1.0
  #function fb(x,y)
  #  rad=sqrt((x-xM)^2+(y-zM)^2);
  #  return 0.0+(rad<1000)*(2.0*cos(0.5*pi*rad/1000)^2);
  #end
  #fb(x,y)=(x>60000 && x<140000)*1.0
  #fb(x,y)=1.0;
  f=Dict(:v=>[fb1,fb2])

  assembMass!(p);
  assembStiff!(p);
  applyStartValues!(p, f);

  v1(x,y)=UMax
  v2(x,y)=0.0
  V=[v1, v2];
  Vf=projectAdvection(p,V,Vfcomp);
  #println("Vf")
  #for i in 1:length(Vf)
  #  println("Vf[$i]: $(Vf[i])")
  #end

  taskRecovery ? pos=[1,3] : pos=[1];
  advectionTypes=Symbol[Vfcomp];
  for i in keys(femType)
      append!(advectionTypes,femType[i][pos]);
  end
  nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, collect(Set(advectionTypes)), size(p.kubWeights,2));
  setEdgeData!(p, :v);

  y=p.solution[0.0];
  Time=0.0;
  for i=1:nIter
    y=splitExplicit(p,gamma,Vfcomp,Vf,nquadPhi,nquadPoints,MISMethod,y,Time,dt,ns);
    Time=Time+dt
    p.solution[Time]=y;
    if mod(i,10)==0
      p2=deepcopy(p)
      #unstructured_vtk(p2, sort(collect(keys(p2.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAdvectionTriangles/"*filename)
    end
    println(Time)
  end

  #Speichern des Endzeitpunktes als vtu-Datei:
  #unstructured_vtk(p, EndTime, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAdvectionTriangles/"*filename)
  #Speichern aller berechneten Zwischenwerte als vtz-Datei:
  unstructured_vtk(p, sort(collect(keys(p.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAdvectionTriangles/"*filename)

  return p
end
