include("modulesSphere.jl")

function testSphere()

    filename = "cubedSphereConstant";

    #order: comp, compHigh, compRec, compDG

    femType=Dict(:rho=>[:DG0, :P1, :DG1, :DG0],
                 :rhoV=>[:RT0, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG0, :P1, :DG1, :DG0],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);

    #higher spaces
    #=
    femType=Dict(:rho=>[:DG1, :P1, :DG1, :DG0],
                 :rhoV=>[:RT1, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG1, :P1, :DG1, :DG0],
                 :p=>[:DG1],
                 :v=>[:RT1],
                 :theta=>[:DG1]);
    =#

    taskRecovery=true;
    advection=true;

    m=generateCubedSphere(36,6300000.0)

    p=femProblem(m, femType,t=:shallow, advection=advection, taskRecovery=taskRecovery);
    #return p;
    gamma=0.5; #upwind
    UMax=0.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS2); #method of time integration

    dt=50.0;
    ns=10;
    EndTime=1000.0;
    nIter=Int64(EndTime/dt);

    #start functions
    Rad=6300000.0
    function frho(xyz::Array{Float64,1})
        #=
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        rd=distCircle(lon,lat,0.5*pi,0.0,Rad)
        R=Rad/3.0;
        rd<=R ? rhoLoc=1000.0/2.0*(1.0+cos(pi*rd/R))  : rhoLoc=0.0; #rhoLoc=1000.0/2.0*(1.0+cos(pi*rd/R))
        return rhoLoc+1.0;
        =#
        return 2.0
    end
    function ftheta(xyz::Array{Float64,1})
        return 1.0;
    end
    fvel(xyz::Array{Float64,1})=[UMax,0.0,0.0];
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

    unstructured_vtk3D(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"0")

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    taskRecovery ? pos=[1,3] : pos=[1];
    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        append!(advectionTypes,femType[i][pos]);
    end
    nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, collect(Set(advectionTypes)), size(p.kubWeights,2));
    setEdgeData!(p, :v)

    MrT=assembMass(p.degFBoundary[femType[:rhoTheta][1]], m, p.kubPoints, p.kubWeights);
    MrV=assembMass(p.degFBoundary[femType[:rhoV][1]], m, p.kubPoints, p.kubWeights);

    y=p.solution[0.0];
    Y=Array{solution,1}(undef,MISMethod.nStage+1);
    FY=Array{solution,1}(undef,MISMethod.nStage);
    SthY=Array{SparseMatrixCSC{Float64,Int64},1}(undef,MISMethod.nStage);
    Time=0.0;

    for i=1:nIter
      @time y=splitExplicit(y,Y,FY,SthY,p,gamma,nquadPhi,nquadPoints,MrT,MrV,MISMethod,Time,dt,ns);
      Time+=dt
      p.solution[Time]=y;
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT);
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)
      p2=deepcopy(p);
      unstructured_vtk3D(p2, Time, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"$i")
      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk3D(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk3D(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    return p;
end
