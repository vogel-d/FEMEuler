include("modulesCE.jl")

function testWarmBubble()
    filename = "warmBubble";

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

    m=generateRectMesh(160,80,:periodic,:constant,-10000.0,10000.0,0.0,10000.0); #(east/west, top/bottom)
    #m=generateRectMesh(80,40,:periodic,:constant,-10000.0,10000.0,0.0,10000.0); #(east/west, top/bottom)

    #adaptGeometry!(m,(0.3,0.3),false); #sin perbutation

    p=femProblem(m, femType,t=:compressible, advection=advection, taskRecovery=taskRecovery);

    gamma=0.5; #upwind
    UMax=0.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS2); #method of time integration

    dt=1.0; #Coarse: 2.0
    #dt=0.5; #Coarse: 1.0
    ns=15;
    EndTime=1000.0;
    nIter=Int64(EndTime/dt);

    #start functions
    xCM=0.0; zCM=2000.0;
    r0=2000.0; th0=300.0; p0=100000.0;
    DeltaTh1=2;
    Grav=9.81;
    Cpd=1004.0; Cvd=717.0; Cpv=1885.0;
    Rd=Cpd-Cvd; Gamma=Cpd/Cvd; kappa=Rd/Cpd;
    function frho(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        pLoc=p0*(1-kappa*Grav*z/(Rd*th0))^(Cpd/Rd);
        Rad=sqrt((x-xCM)^2+(z-zCM)^2);
        ThLoc=th0+(Rad<r0)*(DeltaTh1*cos(0.5*pi*Rad/r0)^2);
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function ftheta(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        rad=sqrt((x-xCM)^2+(z-zCM)^2);
        return th0+(rad<r0)*(DeltaTh1*cos(0.5*pi*rad/r0)^2);
    end
    fv1(xz::Array{Float64,1})=UMax;
    fv2(xz::Array{Float64,1})=0.0;
    fvel=[fv1, fv2];
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    p.boundaryValues[(:theta,:P1)]=300.0*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);

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
      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    return p
end
