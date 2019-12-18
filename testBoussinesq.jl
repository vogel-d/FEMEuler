include("modulesB.jl")

function testBoussinesq()
    filename = "testPeriodic";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    boundaryCondition = (:periodic,:constant); #(east/west, top/bottom)

    pv=femProblem(:quad, 300, 10, femType, boundaryCondition, xr=300000.0, yr=10000.0);


    method=:euler;
    dt=1.0;
    tend=3000.0;

    #solSaves=15.0:15:tend; #determines at which points of time the solution is saved
    solSaves=tend;

    b0=0.01;
    H=10000;
    A=5000;
    xM=0.75*(pv.mesh.geometry.l[1]+pv.mesh.geometry.r[1]);
    #yM=0.5*(pv.mesh.geometry.l[2]+pv.mesh.geometry.r[2]);

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
        println(i);
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(pv, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)

    return pv
end
