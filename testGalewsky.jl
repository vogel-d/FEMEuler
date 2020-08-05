include("modulesSphere.jl")

function testGalewsky()

    filename = "galewskyII";

    stencilOrder=2;
    recoveryOrder=2;

    recoverySpace=Symbol("R$recoveryOrder")
    recoverySpaceVec=Symbol("VecR$(recoveryOrder)S")

    #order: comp, compTest, recoverySpace
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 #:rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 :rhoV=>[:RT0, :VecP1, :VecDG1, :RT0B],
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
    femType=Dict(:rho=>[:DG1, :P1, :DG1, :DG0],
                 :rhoV=>[:RT1, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG1, :P1, :DG1, :DG0],
                 :p=>[:DG1],
                 :v=>[:RT1],
                 :theta=>[:DG1]);
    =#

    taskRecovery=true;
    adv=true;

    m=generateCubedSphere(80,6300000.0,0,:cube1)

    p=femProblem(m, femType,t=:shallow, advection=adv, taskRecovery=taskRecovery,
    stencilOrder=stencilOrder, recoveryOrder=recoveryOrder,);

    gamma=0.5; #upwind
    MISMethod=MIS(:MIS4_4); #method of time integration

    dt=200.0 #200.0 #50.0 #1200.0;
    ns=120;
    EndTime=6.0*86400.0;
    nIter=Int64(EndTime/dt);

    #start functions
    alphaG=1.0/3.0
    betaG=1.0/15.0
    hH=120.0
    H0G=10000.0
    function integrandG(tau,RadEarth)
      uM=80.0
      lat0G=pi/7.0
      lat1G=pi/2.0-lat0G
      eN=exp(-4.0/(lat1G-lat0G)^2.0)
      f=2.0*Omega*sin(tau)
      if tau<=lat0G || tau>=lat1G
        uStart=0.0
      else
        uStart=uM/eN*exp(1.0/((tau-lat0G)*(tau-lat1G)))
      end
      if abs(tau)<0.5*pi
        intG=(RadEarth*f+uStart*tan(tau))*uStart
      else
        intG=0.0
      end
      return intG
    end
    function frho(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        return (Grav*H0G-(simpson(-0.5*pi,lat,r,pi/100.0,integrandG)))/Grav#+hH*cos(lat)*exp(-((lon-pi)/alphaG)^2.0)*exp(-((pi/4.0-lat)/betaG)^2.0)
    end
    function ftheta(xyz::Array{Float64,1})
        return 1.0;
    end
    function fvel(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        uM=80.0
        lat0G=pi/7.0
        lat1G=pi/2.0-lat0G
        eN=exp(-4.0/(lat1G-lat0G)^2.0)
        if lat<=lat0G || lat>=lat1G
          uS=0
        else
          uS=uM/eN*exp(1.0/((lat-lat0G)*(lat-lat1G)))
        end
        return velCa([uS,0.0,0.0],lon,lat)
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"0", printSpherical=true)

    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        push!(advectionTypes,femType[i][1]);
        (taskRecovery && length(femType[i])==4) && push!(advectionTypes,femType[i][3]);
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

    Sth=advection(p,gamma,y,nquadPoints,nquadPhi,MrT)
    for i=1:nIter
      @time y=splitExplicit(y,Y,FY,SthY,p,gamma,Sth,nquadPhi,nquadPoints,MrT,MrV,MISMethod,Time,dt,ns);
      Time+=dt
      p.solution[Time]=y;
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT);
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)
      if mod(i,8)==0
          p2=deepcopy(p);
          unstructured_vtk(p2, Time, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"$i", printSpherical=true)
      end
      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename, printSpherical=true)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)

    return p;
end
p=testGalewsky()
