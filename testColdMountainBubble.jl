include("modulesCE.jl")

const stencilOrder=2;

recoverySpace=Symbol("DGQuad")
recoverySpaceVec=Symbol("VecDGQuad")

@recovery(recoverySpace,recoverySpaceVec)

function testColdMountainBubble()
    filename = "coldMountainBubbleTRQuad";

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 :rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 :rhoTheta=>[:DG0, :DG0, recoverySpace],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);
    #=
    femType=Dict(:rho=>[:DG0, :P1, :DG1, :DG0],
                 :rhoV=>[:RT0, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG0, :P1, :DG1, :DG0],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);
    =#
    #higher spaces
    #=
    femType=Dict(:rho=>[:DG1, :DG1, recoverySpace],
                 :rhoV=>[:RT1, :RT1, recoverySpaceVec],
                 :rhoTheta=>[:DG1, :DG1, recoverySpace],
                 :p=>[:DG1],
                 :v=>[:RT1],
                 :theta=>[:DG1]);
    =#

    taskRecovery=true;
    advection=true;

    m=generateRectMesh(360,64,:periodic,:constant,-18000.0,18000.0,0.0,6400.0); #(east/west, top/bottom)
    #m=generateRectMesh(180,32,:periodic,:constant,-18000.0,18000.0,0.0,6400.0); #(east/west, top/bottom)
    p=femProblem(m, femType, t=:compressible, advection=advection, taskRecovery=taskRecovery,
    stencilOrder=stencilOrder);

    adaptGeometry!(m,1000.0,2000.0); #witch of agnesi with Gall-Chen and Sommerville transformation

    gamma=0.5; #upwind
    UMax=0.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS4_4); #method of time integration

    dt=0.5 #Coarse: 1.0, Fine: 0.5
    ns=6;
    EndTime=900.0;
    nIter=Int64(EndTime/dt)

    #start functions
    xCM=0.0; zCM=3000.0;
    xCR=4000.0; zCR=2000.0;
    r0=2000.0; th0=300.0;
    DeltaTh1=-15;
    function frho(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        pLoc=p0*(1.0-kappa*Grav*z/(Rd*th0))^(Cpd/Rd)
        Rad=sqrt(((x-xCM)/xCR)^2+((z-zCM)/zCR)^2);
        ThStart=th0
        if Rad<1.0
            ThStart+=DeltaTh1*(cos(pi*Rad)+1.0)/2.0*(pLoc/p0)^(-kappa)
        end
        return pLoc/((pLoc/p0)^kappa*Rd*ThStart)
    end
    function ftheta(xz::Array{Float64,1})
        x=xz[1]; z=xz[2];
        pLoc=p0*(1.0-kappa*Grav*z/(Rd*th0))^(Cpd/Rd)
        Rad=sqrt(((x-xCM)/xCR)^2+((z-zCM)/zCR)^2);
        ThStart=th0
        if Rad<1.0
            ThStart+=DeltaTh1*(cos(pi*Rad)+1.0)/2.0*(pLoc/p0)^(-kappa)
        end
        return ThStart
    end
    function fvel(xz::Array{Float64,1})
        return [UMax, 0.0]
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    #p.boundaryValues[(:theta,:P1)]=300*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"0")

    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        push!(advectionTypes,femType[i][1]);
        push!(advectionTypes,femType[i][3]);
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
      Time=Time+dt
      p.solution[Time]=y;
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT);
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)

      if mod(i,50)==0
          p2=deepcopy(p);
          unstructured_vtk(p2, Time, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"$i")
      end

      println(Time)
    end
    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, collect(keys(p.solution)), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    return p
end
