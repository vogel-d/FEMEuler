include("../src/Modules/modulesB.jl")

function testBoussinesq()
    filename = "boussinesqCG1";

    femType=Dict(:p=>[:DG0], :v=>[:VecP1], :b=>[:DG0]);
    #femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    pv=femProblem(m, femType);

    method=:rk4  #:euler;
    dt=1.0; #Fine: 0.5
    tend=3000.0;

    #determines at which points of time the solution is saved
    #solSaves=1;
    solSaves=10;
    nIter=tend/solSaves;

    b0=0.01;
    H=10000;
    A=5000;
    xM=0.5*(pv.mesh.geometry.l[1]+pv.mesh.geometry.r[1]);

    function fb(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        return b0*sin(pi*z/H)/(1+((x-xM)/A)^2);
    end
    f=Dict(:b=>fb);

    assembMass!(pv);
    assembStiff!(pv);
    applyStartValues!(pv,f)

    Fp=pv.massM[pv.femType[:p][1]];
    Fv=pv.massM[pv.femType[:v][1]];
    Fb=pv.massM[pv.femType[:b][1]];

    time=0.0
    for i in 1:solSaves
        solveB!(pv,Fp,Fv,Fb,time,dt,nIter,method);
        time+=nIter
        println(time)
    end
    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(pv, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)

    return pv
end
p=testBoussinesq();
