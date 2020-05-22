include("modulesSphereAdv.jl")

function testSphereAdv()

    filename = "testSphereAdvTRSH";

    stencilOrder=2;
    recoveryOrder=2;

    recoverySpace=Symbol("R$recoveryOrder")
    recoverySpaceVec=Symbol("VecR$recoveryOrder")

    #order: comp, compTest, recoverySpace
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 #:rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 :rhoV=>[:RT0, :VecP1S, :VecDG1S, :RT0B],
                 :rhoTheta=>[:DG0, :DG0, recoverySpace],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);
    #=
    femType=Dict(:rho=>[:DG0, :P1, :DG1, :DG0],
                 :rhoV=>[:RT0, :VecP1S, :VecDG1S, :RT0B],
                 :rhoTheta=>[:DG0, :P1, :DG1, :DG0],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :theta=>[:DG0]);

    #higher spaces

    femType=Dict(:rho=>[:DG1, :P1, :DG1, :DG0],
                 :rhoV=>[:RT1, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG1, :P1, :DG1, :DG0],
                 :p=>[:DG1],
                 :v=>[:RT1],
                 :theta=>[:DG1]);
    =#

    taskRecovery=true;
    adv=true;

    m=generateCubedSphere(15,sqrt(3.0),0,:cube1)

    p=femProblem(m, femType,t=:shallow, advection=adv,
        taskRecovery=taskRecovery, stencilOrder=stencilOrder,
        recoveryOrder=recoveryOrder, g=4);

    gamma=0.5; #upwind
    UMax=5.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MISRK2); #method of time integration

    dt=0.0105;
    ns=15;
    #EndTime=40000.0
    #nIter=Int64(EndTime/dt);
    nIter=200;

    #start functions
    function frho(xyz::Array{Float64,1})
        return 1.0
    end
    function ftheta(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];

        #y+
        lat0=4.0*atan(1.0)
        lon0=2.0*atan(1.0)

        lon,lat,r=cart2sphere(x,y,z);
        d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
        if abs(d)<=0.8 #0.4
            conc=1.0
        else
            conc=0.1
        end
        return conc;
    end

    function fvel(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        #=
        if caseVel==:Spherical
            return velCa([UMax*cos(lat),0.0,0.0],lon,lat)
        elseif caseVel==:BlobC
            lat0=4.0*atan(1.0)
            lon0=1.8*atan(1.0) #1.25*atan(1.0) #1.5*atan(1.0)
            d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
            if abs(d)<=0.4 #0.8 #0.1
                uS=1.0
            else
                uS=0.1
            end
            return [uS,0.0,0.0]
        elseif caseVel==:BlobS
        =#
            lat0=4.0*atan(1.0)
            lon0=1.8*atan(1.0) #1.25*atan(1.0) #1.5*atan(1.0)
            d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
            if abs(d)<=0.4 #0.8 #0.1
                uS=1.0
            else
                uS=0.1
            end
            return velCa([uS,0.0,0.0],lon,lat)
        #end
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    #assembStiff!(p);
    applyStartValues!(p, f);

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"0", printSpherical=true)

    function V(xyz::Array{Float64,1})
      x=xyz[1]; y=xyz[2]; z=xyz[3];
      λ,θ,r=cart2sphere(x,y,z);
      α=0.0 #0.0, 0.05, pi/2-0.05, pi/2
      uS=UMax*(cos(θ)*cos(α)+sin(θ)*cos(λ)*sin(α))
      vS=-UMax*sin(λ)*sin(α)
      return velCa([uS,vS,0.0],λ,θ)
    end

    Vfcomp=:RT0
    Vf=projectAdvection(p,V,Vfcomp);
    #vtk(m,p.degFBoundary[Vfcomp],Vector(Vf),Vfcomp,"testVf", printSpherical=true)

    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        push!(advectionTypes,femType[i][1]);
        (taskRecovery && length(femType[i])==4) && push!(advectionTypes,femType[i][3]);
    end
    nquadPhi, nquadPoints=coordTrans(m, m.normals, advectionTypes, size(p.kubWeights,2));
    setEdgeData!(p, :v)

    MrT=assembMass(p.degFBoundary[femType[:rhoTheta][1]], m, p.kubPoints, p.kubWeights);
    MrV=assembMass(p.degFBoundary[femType[:rhoV][1]], m, p.kubPoints, p.kubWeights);

    y=p.solution[0.0];
    Y=Array{solution,1}(undef,MISMethod.nStage+1);
    FY=Array{solution,1}(undef,MISMethod.nStage);
    SthY=Array{SparseMatrixCSC{Float64,Int64},1}(undef,MISMethod.nStage);
    Time=0.0;
    for i=1:nIter
      ry=advection(p,gamma,y,Vf,Vfcomp,nquadPoints,nquadPhi,MrT,MrV);
      yNeu=y+0.5*dt*ry
      ry=advection(p,gamma,yNeu,Vf,Vfcomp,nquadPoints,nquadPhi,MrT,MrV);
      y+=dt*ry;
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
    #=
    for i=1:nIter
      @time y=splitExplicit(y,Y,FY,SthY,p,Vfcomp,Vf,gamma,nquadPhi,nquadPoints,MrT,MrV,MISMethod,Time,dt,ns);
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
    =#
    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename, printSpherical=true)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(0.0:200.0:EndTime)), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)

    return p;
end
