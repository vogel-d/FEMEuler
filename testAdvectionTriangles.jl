include("modulesBA.jl")

function testAdvectionTri()
  filename = "testBubbleY"

  #order: comp, compHigh, compRec, compDG
  femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG0, :P1, :DG1, :DG0]);
  #femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:DG1, :P1, :DG1, :DG0]);
  #femType=Dict(:p=>[:DG0, :P1, :DG1, :DG0], :v=>[:RT0, :VecP1, :VecDG1, :RT0B], :b=>[:P1, :P1, :DG1, :DG1]);
  #femType=Dict(:p=>[:DG1, :P1, :DG1, :DG0], :v=>[:RT1, :VecP1, :VecDG1, :RT0B], :b=>[:DG1, :P1, :DG1, :DG0]);
  Vfcomp=:RT0
  #Vfcomp=:RT1

  taskRecovery=true;

  #m=generateTriMesh(300,20,:periodic,:periodic,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
  m=generateTriMesh(160,80,:periodic,:periodic,-10000.0,10000.0,0.0,10000.0); #(east/west, top/bottom)
  #m=generateTriMesh(3,3,:periodic,:periodic,0.0,6.0,0.0,6.0); #(east/west, top/bottom)
  #adaptGeometry!(m,(0.3,0.3),false); #sin perbutation
  p=femProblem(m, femType, t=:boussinesq, taskRecovery=taskRecovery);

  gamma=0.5;
  UMax=20.0;
  #UMax=1.0; #UMax determines the advection in x direction
  #UMax=0.0;
  #MISMethod=MIS(:MIS_Euler);
  MISMethod=MIS(:MIS4_4);

  #dt=10.0;
  dt=1.0;
  #dt=0.5;
  ns=19;
  #EndTime=1*dt;
  EndTime=1000.0;
  nIter=Int64(EndTime/dt);

  #start function
  xR=m.geometry.r[1]; xL=m.geometry.l[1]; zR=m.geometry.r[2]; zL=m.geometry.l[2]
  b0=0.01; H=10000; A=5000; r0=2000.0;
  xM=0.5*(xL+xR); zM=0.5*(zL+zR);
  #function fb1(xz::Array{Float64,1})
  #      x=xz[1]; z=xz[2];
  #      return b0*sin(pi*z/H)/(1+((x-xM)/A)^2);
  #end
  #fb2(x,y)=1.0
  function fb1(xz::Array{Float64,1})
    x=xz[1]; z=xz[2];
    rad=sqrt((x-xM)^2+(z-zM)^2);
    return 0.0+(rad<r0)*(2.0*cos(0.5*pi*rad/r0)^2);
  end
  function fb2(xz::Array{Float64,1})
    x=xz[1]; z=xz[2];
    return (x>1 && x<2 && z>1 && z<2)*1.0
  end
  function fv1(xz::Array{Float64,1})
    return 1
  end
  function fv2(xz::Array{Float64,1})
    return 0
  end
  #fb(x,y)=1.0;
  #f=Dict(:b=>fb1)
  f=Dict(:v=>[fv2,fb1])

  assembMass!(p);
  assembStiff!(p);
  applyStartValues!(p, f);

  v1(xz::Array{Float64,1})=UMax;
  v2(xz::Array{Float64,1})=0.0;
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

  gamma=0.5;

  #time-stepping
  RKadvection!(p,gamma,Vfcomp,Vf,nquadPhi,nquadPoints,f,Time,dt,EndTime);



#=
  for i=1:nIter
    #y=splitExplicit(p,gamma,Vfcomp,Vf,nquadPhi,nquadPoints,MISMethod,f,Time,dt);
    y=RKadvection!(p,gamma,Vfcomp,Vf,nquadPhi,nquadPoints,f,Time,dt);
    Time=Time+dt
    p.solution[Time]=y;
    if mod(i,10)==0
      #p2=deepcopy(p)
      #unstructured_vtk(p2, sort(collect(keys(p2.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAdvectionTriangles/"*filename)
    end
    println(Time)
  end
=#

  #Speichern des Endzeitpunktes als vtu-Datei:
  #unstructured_vtk(p, EndTime, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAdvectionTriangles/"*filename)
  #Speichern aller berechneten Zwischenwerte als vtz-Datei:
  unstructured_vtk(p, sort(collect(keys(p.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAdvectionTriangles/"*filename)

  return p
end
