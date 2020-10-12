include("modulesBCompound.jl")

function testCompoundBoussinesq()
    filename = "fail";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:DG0]);
    #femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    #m=generateHexMesh(0.0,300000.0,0.0,10000.0,10,:periodic,:constant,meshType=3); #(east/west, top/bottom)
    #m=generateHexMesh(0.0,300000.0,0.0,10000.0,10,:periodic,:constant,meshType=4); #(east/west, top/bottom)
    m=generateHexMesh(100000.0,200000.0,0.0,10000.0,10,:periodic,:constant,meshType=4); #(east/west, top/bottom)
    #m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    #m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0,meshType=3); #(east/west, top/bottom)
    #pv=femProblem(m, femType, compoundMethod=:HexToTris);
    #pv=femProblem(m, femType, compoundMethod=:RectToTris);
    pv=femProblem(m, femType, compoundMethod=:HexToKites);
    #pv=femProblem(m, femType, compoundMethod=:RectToKites);

    #adaptGeometry!(pv.mesh,100.0);

    assemblePhiPre!(pv);

    method=:euler;
    dt=1.0; #Fine: 0.5
    tend=dt;
    #tend=400.0;
    #tend=3000.0;

    #determines at how many points of time the solution is saved
    #solSaves=100;
    solSaves=1;
    nIter=max(tend/solSaves,1);

    b0=0.01;
    #H=10000;
    H=pv.mesh.geometry.r[2];
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
    compound_unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBoussinesq/"*"compound_"*filename)
    #compound_unstructured_vtk(pv, [0.0,500.0,1000.0,3000.0], [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBoussinesq/"*"compound_"*filename)
    #split_unstructured_vtk(pv, sort(collect(keys(pv.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testCompoundBoussinesq/"*"split_"*filename)

    return pv
end
