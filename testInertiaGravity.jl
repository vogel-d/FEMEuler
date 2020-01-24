include("modulesCEd.jl")

function testInertiaGravity()
    filename = "gravityWavesNoAdv";

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :P1, :DG1, :DG0],
                 :rhoV=>[:RT0, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG0, :P1, :DG1, :DG0],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0],
                 :pBar=>[:DG0],
                 :rhoBar=>[:DG0],
                 :thBar=>[:DG0]);

    taskRecovery=false;
    advection=true;

    m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    p=femProblem(m, femType, t=:compressible, advection=advection, taskRecovery=taskRecovery);

    gamma=0.5; #upwind
    UMax=20.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS4_4); #method of time integration

    dt=10.0;
    ns=10;
    EndTime=3000.0;
    nIter=div(EndTime,dt);

    #start functions
    xCM=0.0; zCM=2000.0;
    r0=2000.0; th0=300.0; p0=100000.0;
    DeltaTh1=.01;
    Grav=9.81;
    Cpd=1004.0; Cvd=717.0; Cpv=1885.0;
    N=1.e-2;
    Rd=Cpd-Cvd; Gamma=Cpd/Cvd; kappa=Rd/Cpd;
    H=10000;
    a=5000;
    xC=150000;
    function frho(x::AbstractFloat,z::AbstractFloat)
        S=N*N/Grav
        ThLoc=th0*exp(S*z)+DeltaTh1*sin(pi*z/H)/(1.0+((x-xC)/a)^2);
        pLoc=p0*(1.0-Grav/(Cpd*th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd)
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function frhoBar(x::AbstractFloat,z::AbstractFloat)
        S=N*N/Grav
        ThLoc=th0*exp(S*z)
        pLoc=p0*(1.0-Grav/(Cpd*th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd)
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function fpBar(x::AbstractFloat,z::AbstractFloat)
        S=N*N/Grav
        return p0*(1.0-Grav/(Cpd*th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd)
    end
    function fthBar(x::AbstractFloat,z::AbstractFloat)
        S=N*N/Grav
        return th0*exp(S*z)
    end
    function ftheta(x::AbstractFloat,z::AbstractFloat)
        S=N*N/Grav
        return th0*exp(S*z)+DeltaTh1*sin(pi*z/H)/(1.0+((x-xC)/a)^2);
    end
    fv1(x, y)=UMax;
    fv2(x, y)=0;
    fvel=[fv1, fv2];
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel,:rhoBar=>frhoBar,:pBar=>fpBar,:thBar=>fthBar);

    assembMass!(p);
    assembStiff!(p);
    p.boundaryValues[(:theta,:P1)]=300*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);

    p.solution[0.0].theta-=p.diagnostic.thBar;
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
    SthY=Array{SparseMatrixCSC{AbstractFloat,Int},1}(undef,MISMethod.nStage);
    Time=0.0;
    for i=1:nIter
      @time y=splitExplicit(y,Y,FY,SthY,p,gamma,nquadPhi,nquadPoints,MrT,MrV,MISMethod,Time,dt,ns);
      Time+=dt
      p.solution[Time]=y;
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT)#-p.diagnostic.thBar;
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)

      if mod(i,50)==0
        p2=deepcopy(p);
        unstructured_vtk(p2, Time, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
      end

      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    #unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"Start")
    return p
end
