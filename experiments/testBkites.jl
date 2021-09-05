include("../src/Modules/modulesB.jl")
include("../src/Transformations/coordTrans.jl")

function testBoussinesq()
    filename = "bachelor_boussinesq_kite";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:DG0]);
    #femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    #m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    #pv=femProblem(m, femType);

    mHex=generateHexMesh(0.0,300000.0,0.0,10000.0,10,:periodic,:constant,meshType=4); #(east/west, top/bottom)
    p=femProblem(mHex, femType, compoundMethod=:HexToKites);
    m=splitCompoundMeshKites(p);
    pv=femProblem(m, femType);



    method=:rk4  #:euler;
    dt=0.25; #Fine: 0.5
    tend=3000.0;

    #determines at how many points of time the solution is saved
    solSaves=30;
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
    #unstructured_vtk(pv, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)
    unstructured_vtk(pv, [0.0,500.0,1000.0,3000.0], [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesq/"*filename)

    return pv
end
