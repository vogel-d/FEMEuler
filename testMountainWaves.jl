include("modulesCE.jl")

function testMountainWaves()
    filename = "mountainWaves";

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :P1, :DG1, :DG0],
                 :rhoV=>[:RT0, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG0, :P1, :DG1, :DG0],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);

    taskRecovery=true;
    advection=true;

    #p=femProblem(:quad, 200, 156, femType, t=:compressible, advection=advection, taskRecovery=taskRecovery,
    #            xl=-20000.0,xr=20000.0, yr=15600.0);
    p=femProblem(:quad, 200, 90, femType, t=:compressible, advection=advection, taskRecovery=taskRecovery,
                xl=-20000.0,xr=20000.0, yr=9000.0);
    println("Problem erstellt")
    adaptGeometry!(p,400.0,1000.0); #witch of agnesi with Gall-Chen and Sommerville transformation
    println("Mesh verzogen")
    boundaryCondition = (:periodic, :constant); #(top/bottom, east/west)

    gamma=0.5; #upwind
    UMax=10.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS2); #method of time integration

    dt=10.0;
    ns=20;
    EndTime=2160.0;
    nIter=Int64(EndTime/dt);
    nIter=5;

    #start functions
    th0=300.0; p0=100000.0;
    Grav=9.81; N=0.01
    Cpd=1004.0; Cvd=717.0; Cpv=1885.0;
    Rd=Cpd-Cvd; Gamma=Cpd/Cvd; kappa=Rd/Cpd;
    function frho(x::Float64,z::Float64)
        s=N*N/Grav
        ThLoc=th0*exp(z*s)
        pLoc=p0*(1-Grav/(Cpd*th0*s)*(1-exp(-s*z)))^(Cpd/Rd)
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function ftheta(x::Float64,z::Float64)
        return th0*exp(z*N*N/Grav)
    end

    fv1(x::Float64, y::Float64)=UMax;
    fv2(x::Float64, y::Float64)=0.0;
    fvel=[fv1, fv2];
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembFEM!(p, boundaryCondition);
    println("Matrizen berechnet")
    p.boundaryValues[(:theta,:P1)]=300*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);
    println("Startwerte proijziert")

    taskRecovery ? pos=[1,3] : pos=[1];
    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        append!(advectionTypes,femType[i][pos]);
    end
    nquadPhi, nquadPoints=coordTrans(p.mesh.meshType, p.mesh.normals, collect(Set(advectionTypes)), size(p.kubWeights,2));
    setEdgeData!(p, :v)

    MrT=assembMass(p.degFBoundary[femType[:rhoTheta][1]], p.mesh, p.kubPoints, p.kubWeights);
    MrV=assembMass(p.degFBoundary[femType[:rhoV][1]], p.mesh, p.kubPoints, p.kubWeights);

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
      #=
      if mod(i,50)==0
        p2=deepcopy(p);
        unstructured_vtk(p2, sort(collect(keys(p2.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
      end
      =#
      println(Time)
    end
    correctVelocity!(p);
    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    return p

    return m;
end

function adapt!(m::mesh,hm, a)
    H=m.geometry.r[2];
    #hm=1.5
    #a=5
    h(x)=(hm*a^2)/(x^2+a^2);
    coord=m.geometry.coordinates;
    for i in 1:(m.topology.n[1]+1)
        coord[:,i]=[coord[1,i], h(coord[1,i])];
    end
    for i in (m.topology.n[1]+2):size(coord,2)
        z0=h(coord[1,i])
        coord[:,i]=[coord[1,i], H*(coord[2,i]+z0)/(H+z0)];
    end

    return nothing;
end
