include("../src/Modules/modulesCE.jl")

const stencilOrder=2;

recoverySpace=Symbol("DGQuad")
recoverySpaceVec=Symbol("VecDGQuad")

@recovery(recoverySpace,recoverySpaceVec)

function testInertiaGravity2()
    filename = "gravityWavesNoAdvNoBackground";

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 :rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 :rhoTheta=>[:DG0, :DG0, recoverySpace],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);

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

    function ftheta(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        S=N*N/Grav
        return th0*exp(S*z)+DeltaTh1*sin(pi*z/H)/(1.0+((x-xC)/a)^2);
    end
    function fvel(xz::Array{Float64,1})
        return [UMax, 0.0]
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

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
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT);
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)

      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtk-Datei:
    unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    #unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"Start")
    return p
end
p=testInertiaGravity2();

#Startbedingungen ohne Background
#=
#start functions
xW=500.0; zW=300.0; xC=560.0; zC=640.0;
rW0=150.0; rC0=0.0;
th0=300.0; p0=100000.0;
DeltaThW=0.5; DeltaThC=-0.15;
sW=50; sC=50;
Grav=9.81;
Cpd=1004.0; Cvd=717.0; Cpv=1885.0;
Rd=Cpd-Cvd; Gamma=Cpd/Cvd; kappa=Rd/Cpd;
function frho(x::Float64,z::Float64)
    pLoc=p0*(1-kappa*Grav*z/(Rd*th0))^(Cpd/Rd);
    #Rad=sqrt((x-xCM)^2+(z-zCM)^2);
    #ThLoc=th0+(Rad>r0)*(DeltaTh1*exp(-(Rad-r0)^2/s^2));
    radW=sqrt((x-xW)^2+(z-zW)^2);
    radC=sqrt((x-xC)^2+(z-zC)^2);
    ThLoc=th0;
    if radW>rW0
        ThLoc+=DeltaThW*exp(-(radW-rW0)^2/sW^2)
    else
        ThLoc+=DeltaThW
    end
    if radW>rC0
        ThLoc+=DeltaThC*exp(-(radC-rC0)^2/sC^2)
    else
        ThLoc+=DeltaThC
    end
    return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
end
function ftheta(x::Float64,z::Float64)
    radW=sqrt((x-xW)^2+(z-zW)^2);
    radC=sqrt((x-xC)^2+(z-zC)^2);
    th=th0;
    if radW>rW0
        th+=DeltaThW*exp(-(radW-rW0)^2/sW^2)
    else
        th+=DeltaThW
    end
    if radW>rC0
        th+=DeltaThC*exp(-(radC-rC0)^2/sC^2)
    else
        th+=DeltaThC
    end
    return th;
    #return th0+(radW>r0)*(DeltaTh1*exp(-(rad-r0)^2/s^2))+(radC>r0)*(DeltaTh1*exp(-(rad-r0)^2/s^2));
end
=#
