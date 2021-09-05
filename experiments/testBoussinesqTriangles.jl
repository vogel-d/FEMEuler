include("../src/Modules/modulesB.jl")

function testBoussinesqTri()
    filename = "testjez";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:P1]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    m=generateTriMesh(300,20,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    #m=generateTriMesh(2,2,:periodic,:constant,0.0,1.0,0.0,1.0); #(east/west, top/bottom)
    p=femProblem(m, femType, g=5);


    method=:euler;
    dt=0.5;
    tend=3000.0;

    solSaves=100; #determines every which point of time the solution is saved
    nIter=tend/solSaves;
    #solSaves=tend;

    b0=0.01;
    H=10000;
    A=5000;
    xM=0.5*(p.mesh.geometry.l[1]+p.mesh.geometry.r[1]);
    #yM=0.5*(p.mesh.geometry.l[2]+p.mesh.geometry.r[2]);

    function fb(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        return b0*sin(pi*z/H)/(1+((x-xM)/A)^2);
    end
    f=Dict(:b=>fb);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p,f)

    Fp=p.massM[p.femType[:p][1]];
    Fv=p.massM[p.femType[:v][1]];
    Fb=p.massM[p.femType[:b][1]];

    #return p;
    time=0.0;
    for i in 1:solSaves
        solveB!(p,Fp,Fv,Fb,time,dt,nIter,method);
        time+=nIter;
    end
    println(time)
    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqTriangles/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    unstructured_vtk(p, sort(collect(keys(p.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqTriangles/"*filename)

    return p
end
