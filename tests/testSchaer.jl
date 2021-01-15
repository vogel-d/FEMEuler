include("modulesCE.jl")

const stencilOrder=2;

recoverySpace=Symbol("DGQuad")
recoverySpaceVec=Symbol("VecDGQuad")

@recovery(recoverySpace,recoverySpaceVec)

function testSchaer()
    filename = "mountainSchaerTRQuad";

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 :rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 :rhoTheta=>[:DG0, :DG0, recoverySpace],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);

    taskRecovery=true;
    advection=true;

    m=generateRectMesh(400,65,:periodic,:constant,-100000.0,100000.0,0.0,19500.0); #(east/west, top/bottom)

    adaptGeometry!(m,250.0,5000.0,4000.0); #Schaer mountain with SLEVE transformation

    p=femProblem(m, femType, t=:compressible, advection=advection, taskRecovery=taskRecovery,
    stencilOrder=stencilOrder);

    gamma=0.5; #upwind
    UMax=10.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS2); #method of time integration

    dt=3.0; #Evtl. kleiner
    ns=20;
    EndTime=2160.0;
    nIter=Int64(EndTime/dt);

    #start functions
    th0=288.0;
    function frho(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        s=N*N/Grav
        ThLoc=th0*exp(z*s)
        pLoc=p0*(1-Grav/(Cpd*th0*s)*(1-exp(-s*z)))^(Cpd/Rd)
        return pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
    end
    function ftheta(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        return th0*exp(z*N*N/Grav)
    end
    function fvel(xz::Array{Float64,1})
        return [UMax, 0.0]
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    println("Matrizen berechnet")
    #p.boundaryValues[(:theta,:P1)]=300*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);
    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);
    println("Startwerte proijziert")

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"0")

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
      y=splitExplicit(y,Y,FY,SthY,p,gamma,nquadPhi,nquadPoints,MrT,MrV,MISMethod,Time,dt,ns);
      Time+=dt
      p.solution[Time]=y;
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT);
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)
      if mod(i,50)==0
        p2=deepcopy(p);
        unstructured_vtk(p2, maximum(collect(keys(p2.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"$i")
      end
      println(Time)
    end
    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    return p
end
