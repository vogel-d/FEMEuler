include("modulesB.jl")

function testBoussinesq()
    filename = "boussinesq";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:DG0]);
    #femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    pv=femProblem(m, femType);

    method=:rk4  #:euler;
    dt=1.0;
    tend=3000.0;

    #determines at which points of time the solution is saved
    #solSaves=15.0:15:tend;
    solSaves=tend;

    b0=0.01;
    H=10000;
    A=5000;
    xM=0.5*(pv.mesh.geometry.l[1]+pv.mesh.geometry.r[1]);

    fb(x,y)=b0*sin(pi*y/H)/(1+((x-xM)/A)^2);
    f=Dict(:b=>fb);

    assembMass!(pv);
    assembStiff!(pv);
    applyStartValues!(pv,f)

    Fp=pv.massM[pv.femType[:p][1]];
    Fv=pv.massM[pv.femType[:v][1]];
    Fb=pv.massM[pv.femType[:b][1]];

    for i in collect(solSaves)
        solveB!(pv,Fp,Fv,Fb,dt,i,method);
    end
    println(tend)
    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(pv, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)

    return pv
end
