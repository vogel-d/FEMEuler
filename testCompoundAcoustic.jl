include("modulesBcompound.jl")
include("solveAcoustic.jl")

function testCompoundAcoustic()
    filename = "testAcousticbig";

    femType=Dict(:p=>[:DG0], :v=>[:RT0], :b=>[:DG0]);
    #femType=Dict(:p=>[:DG1], :v=>[:RT1], :b=>[:DG1]);

    m=generateHexMesh(0.0,1.0,0.0,1.0,3,:periodic,:constant); #(east/west, top/bottom)
    #m=generateRectMesh(3,3,:periodic,:periodic,0.0,2.0,0.0,2.0); #(east/west, top/bottom)
    p=femProblem(m, femType, compoundMethod=:HexToKites);
    #p=femProblem(m, femType);

    method=:euler;
    dt=0.5;
    tend=50.0;

    solSaves=Int(tend/dt); #determines at which points of time the solution is saved
    nIter=tend/solSaves;

    b0=0.01;
    H=10000;
    A=5000;
    xM=0.5*(p.mesh.geometry.l[1]+p.mesh.geometry.r[1]);
    #yM=0.5*(p.mesh.geometry.l[2]+p.mesh.geometry.r[2]);
    xCM=50000.0; zCM=50000.0;
    r0=10000.0; th0=300.0; p0=100000.0;
    DeltaTh1=2;
    Grav=9.81;
    Cpd=1004.0; Cvd=717.0; Cpv=1885.0;
    Rd=Cpd-Cvd; Gamma=Cpd/Cvd; kappa=Rd/Cpd;

    function fp(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        rad=sqrt((x-xCM)^2+(z-zCM)^2);
        return 0+(rad<r0)*(DeltaTh1*cos(0.5*pi*rad/r0)^2)
    end

#=
    function fp(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        if x>2 && x<3 && z>2 && z<3
            result=2.0
        else
            result= 1.0
        end

        return result
    end
=#
    f=Dict(:p=>fp);

    assembMassCompound!(p);
    error()
    assembStiff!(p);
    applyStartValues!(p,f)

    Fp=p.massM[p.femType[:p][1]];
    Fv=p.massM[p.femType[:v][1]];
    Fb=p.massM[p.femType[:b][1]];

    #return p;
    time=0.0;
    for i in 1:solSaves
        solveAcoustic!(p,Fp,Fv,Fb,time,dt,nIter,method);
        time+=nIter;
        #println(time);
    end
    println(time)
    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, tend, [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testBoussinesqTriangles/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    unstructured_vtk(p, sort(collect(keys(p.solution))), [:p, :b, :v], ["Pressure", "Buoyancy", "Velocity"], "testAcousticTriangles/"*filename)

    return p
end
