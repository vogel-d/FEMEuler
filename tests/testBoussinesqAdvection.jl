include("modulesBA.jl")
#include("advectionStiffN.jl")

const stencilOrder=1;

const recoverySpace=Symbol("DGLin")
const recoverySpaceVec=Symbol("VecDGLin")

@recovery(recoverySpace,recoverySpaceVec)

function testBoussinesqAdvection()
  filename = "boussinesqAdvectionTR"

  #femType=Dict(:p=>[:DG0, :DG0, recoverySpace], :v=>[:RT0, :RT0, recoverySpaceVec], :b=>[:DG0, :DG0, recoverySpace]);

  #order: comp, compHigh, compRec, compDG
  femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG0, :P1, :DG1, :DG0]);
  #femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:P1, :P1, :DG1, :DG1]);
  #femType=Dict(:p=>[:DG1, :P1, :DG1, :DG0], :v=>[:RT1, :VecP1, :VecDG1, :RT0B], :b=>[:DG1, :P1, :DG1, :DG0]);
  Vfcomp=:RT0
  #Vfcomp=:RT1

  taskRecovery=true;

  m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
  #m=generateRectMesh(3,3,:periodic,:periodic); #(east/west, top/bottom)
  #adaptGeometry!(m,(0.3,0.3),false); #sin perbutation
  p=femProblem(m, femType, taskRecovery=taskRecovery, stencilOrder=stencilOrder);

  gamma=0.5;
  UMax=20.0; #UMax determines the advection in x direction
  #UMax=0.0;
  MISMethod=MIS(:MIS4_4);

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
  function fb2(xz::Array{Float64,1})
    x=xz[1]; z=xz[2];
    return (x>1 && x<2 && z>1 && z<2)*1.0
  end
  f=Dict(:b=>fb)

  assembMass!(p);
  assembStiff!(p);
  applyStartValues!(p, f);

  unstructured_vtk(p, 0.0, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqAdvection/"*filename*"0")

  V(xz::Array{Float64,1})=[UMax, 0.0];
  Vf=projectAdvection(p,V,Vfcomp);

  advectionTypes=Symbol[Vfcomp];
  for i in keys(femType)
      push!(advectionTypes,femType[i][1]);
      taskRecovery && push!(advectionTypes,femType[i][3]);
  end
  nquadPhi, nquadPoints=coordTrans(m, m.normals, advectionTypes, size(p.kubWeights,2));
  setEdgeData!(p, :v)
  recoveryMatrix!(p)

  y=p.solution[0.0];
  Time=0.0;
  for i=1:nIter
    y=splitExplicit(p,gamma,Vfcomp,Vf,nquadPhi,nquadPoints,MISMethod,y,Time,dt,ns);
    Time=Time+dt
    p.solution[Time]=y;
    #mod(i,5)==0 && unstructured_vtk(p, Time, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqAdvection/"*filename*"$i")
    println(Time)
  end

  #Speichern des Endzeitpunktes als vtu-Datei:
  unstructured_vtk(p, EndTime, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqAdvection/"*filename)
  #Speichern aller berechneten Zwischenwerte als vtk-Datei:
  #unstructured_vtk(p, sort(collect(keys(p.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqAdvection/"*filename)

  return p
end
