include("modulesBAcompound.jl")

function testCompoundBA()
  filename = "testiiiUMax0"

  #order: comp, compHigh, compRec, compDG
  femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG0, :P1, :DG1, :DG0]);
  #femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:P1, :P1, :DG1, :DG1]);
  #femType=Dict(:p=>[:DG1, :P1, :DG1, :DG0], :v=>[:RT1, :VecP1, :VecDG1, :RT0B], :b=>[:DG1, :P1, :DG1, :DG0]);
  Vfcomp=:RT0
  #Vfcomp=:RT1

  taskRecovery=false;

  #m=generateHexMesh(0.0,300000.0,0.0,10000.0,10,:periodic,:constant,meshType=3); #(east/west, top/bottom)
  m=generateHexMesh(0.0,300000.0,0.0,10000.0,10,:periodic,:constant,meshType=4); #(east/west, top/bottom)
  #m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
  #m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0,meshType=3); #(east/west, top/bottom)

  #p=femProblem(m, femType, compoundMethod=:HexToTris);
  p=femProblem(m, femType, compoundMethod=:HexToKites);
  #p=femProblem(m, femType, compoundMethod=:RectToKites);
  #p=femProblem(m, femType, compoundMethod=:RectToTris);

  assemblePhiPre!(p);

  gamma=0.5;
  #gamma=0.0;
  UMax=20.0; #UMax determines the advection in x direction
  #UMax=0.0;
  MISMethod=MIS(:MIS4_4);

  #dt=1.0;
  dt=20.0;
  ns=19;
  EndTime=3000.0;
  nIter=Int64(EndTime/dt);

  #start function
  xR=m.geometry.r[1]; xL=m.geometry.l[1]; yR=m.geometry.r[2]; yL=m.geometry.l[2]
  b0=0.01; H=10000; A=5000;
  xM=0.5*(xL+xR);
  function fb(xz::Array{Float64,1})
    x=xz[1]; z=xz[2];
    return b0*sin(pi*z/H)/(1+((x-xM)/A)^2);
  end
  f=Dict(:b=>fb)

  assembMassCompound!(p);
  assembStiffCompound!(p);
  applyStartValuesCompound!(p, f);

  v1(xz)=UMax;
  v2(xz)=0;
  V=[v1, v2];
  Vf=projectAdvectionCompound(p,V,Vfcomp);

  taskRecovery ? pos=[1,3] : pos=[1];
  advectionTypes=Symbol[Vfcomp];
  for i in keys(femType)
      append!(advectionTypes,femType[i][pos]);
  end
  nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, collect(Set(advectionTypes)), size(p.kubWeights,2));
  setCompoundEdgeData!(p, :v);

  y=p.solution[0.0];
  Time=0.0;
  for i=1:nIter
    @time y=splitExplicit(p,gamma,Vfcomp,Vf,nquadPhi,nquadPoints,MISMethod,y,Time,dt,ns);
    Time=Time+dt
    p.solution[Time]=y;
    if mod(i,10)==0
      p2=deepcopy(p)
      compound_unstructured_vtk(p2, sort(collect(keys(p2.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBA/"*filename)
    end
    println(Time)
  end

  #Speichern des Endzeitpunktes als vtu-Datei:
  #compound_unstructured_vtk(p, EndTime, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBA/"*filename)
  #Speichern aller berechneten Zwischenwerte als vtz-Datei:
  compound_unstructured_vtk(p, sort(collect(keys(p.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBA/"*filename)

  return p
end
