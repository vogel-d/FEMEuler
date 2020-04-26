include("modulesBCompound.jl")

function testCompoundBoussinesq()
    filename = "boussinesq";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:DG0]);
    #femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    m=generateHexMesh(0.0,300000.0,0.0,10000.0,10,:periodic,:constant,meshType=3); #(east/west, top/bottom)
    pv=femProblem(m, femType, compoundMethod=:HexToTris);

    method=:rk4  #:euler;
    dt=1.0; #Fine: 0.5
    tend=3000.0;

    #determines at how many points of time the solution is saved
    solSaves=100;
    #solSaves=15;
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

    assembMassCompound!(pv);
    assembStiffCompound!(pv);
    applyStartValuesCompound!(pv,f)

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
    #unstructured_vtk(pv, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBoussinesq/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBoussinesq/"*filename)

    return pv
end
