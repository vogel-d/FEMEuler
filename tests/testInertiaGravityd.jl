include("modulesCEd.jl")

const stencilOrder=2;

recoverySpace=Symbol("DGQuad")
recoverySpaceVec=Symbol("VecDGQuad")

@recovery(recoverySpace,recoverySpaceVec)

function testInertiaGravity()
    filename = "gravityWavesNoAdv";

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 :rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 :rhoTheta=>[:DG0, :DG0, recoverySpace],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0],
                 :pBar=>[:DG0],
                 :rhoBar=>[:DG0],
                 :thBar=>[:DG0]);

    taskRecovery=false;
    advection=true;

    m=generateRectMesh(300,10,:periodic,:constant,0.0,300000.0,0.0,10000.0); #(east/west, top/bottom)
    p=femProblem(m, femType, t=:compressible, advection=advection, taskRecovery=taskRecovery,
    stencilOrder=stencilOrder);

    gamma=0.5; #upwind
    UMax=20.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS4_4); #method of time integration

    dt=10.0;
    ns=10;
    EndTime=3000.0;
    nIter=div(EndTime,dt);

    #start functions
    xCM=0.0; zCM=2000.0;
    r0=2000.0; th0=300.0;
    DeltaTh1=.01;
    H=10000;
    a=5000;
    xC=150000;
    function frho(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        S=N*N/Grav
        ThLoc=th0*exp(S*z)+DeltaTh1*sin(pi*z/H)/(1.0+((x-xC)/a)^2);
        pLoc=p0*(1.0-Grav/(Cpd*th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd)
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function frhoBar(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        S=N*N/Grav
        ThLoc=th0*exp(S*z)
        pLoc=p0*(1.0-Grav/(Cpd*th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd)
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function fpBar(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        S=N*N/Grav
        return p0*(1.0-Grav/(Cpd*th0*S)*(1.0-exp(-S*z)))^(Cpd/Rd)
    end
    function fthBar(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        S=N*N/Grav
        return th0*exp(S*z)
    end
    function ftheta(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        S=N*N/Grav
        return th0*exp(S*z)+DeltaTh1*sin(pi*z/H)/(1.0+((x-xC)/a)^2);
    end
    function fvel(xz::Array{Float64,1})
        return [UMax, 0.0]
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel,:rhoBar=>frhoBar,:pBar=>fpBar,:thBar=>fthBar);

    assembMass!(p);
    assembStiff!(p);
    #p.boundaryValues[(:theta,:P1)]=300*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);

    p.solution[0.0].theta-=p.diagnostic.thBar;
    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        push!(advectionTypes,femType[i][1]);
        taskRecovery && push!(advectionTypes,femType[i][3]);
    end
    nquadPhi, nquadPoints=coordTrans(m, m.normals, advectionTypes, size(p.kubWeights,2));
    setEdgeData!(p, :v)
    recoveryMatrix!(p)

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
