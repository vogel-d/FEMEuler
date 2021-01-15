include("modulesSphereAdv.jl")

function testSphereAdvS()

    filename = "testSphereAdvS";

    case=:zPanel

    stencilOrder=1;

    recoverySpace=Symbol("R$recoveryOrder")
    recoverySpaceVec=Symbol("VecR$(recoveryOrder)S")

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

    #higher spaces

    femType=Dict(:rho=>[:DG1, :P1, :DG1, :DG0],
                 :rhoV=>[:RT1, :VecP1, :VecDG1, :RT0B],
                 :rhoTheta=>[:DG1, :P1, :DG1, :DG0],
                 :p=>[:DG1],
                 :v=>[:RT1],
                 :theta=>[:DG1]);
    =#

    taskRecovery=true;
    advection=true;

    #m=generateCubedSphere(20,6300000.0)
    #m=generateCubedSphere(5,6300000.0)
    m=generateCubedSphere(3,1.0)

    p=femProblem(m, femType,t=:shallow, advection=advection, taskRecovery=taskRecovery, stencilOrder=stencilOrder);
    #return p;
    gamma=0.5; #upwind
    UMax=10.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS2); #method of time integration

    dt=50.0;
    ns=10;
    EndTime=50.0
    nIter=Int64(EndTime/dt);

    #start functions
    Rad=6300000.0
    function frho(xyz::Array{Float64,1})
        return 1.0
    end
    function ftheta(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lat0=4.0*atan(1.0)
        lon0=2.0*atan(1.0) #1.25*atan(1.0) #1.5*atan(1.0)

        #=
        #y-
        lat0=12.0*atan(1.0)
        lon0=6.0*atan(1.0)
        #x+
        lat0=8.0*atan(1.0)
        lon0=4.0*atan(1.0)
        #x-
        lat0=16.0*atan(1.0)
        lon0=8.0*atan(1.0)
        #z+
        lat0=6.0*atan(1.0)
        lon0=3.0*atan(1.0)
        #z-
        lat0=2.0*atan(1.0)
        lon0=0.0*atan(1.0)
        =#

        #r=sqrt(x*x+y*y+z*z)
        #lat=asin(z/r)
        #lon=atan(x,y)
        lon,lat,r=cart2sphere(x,y,z);
        d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
        if abs(d)<=0.1 #0.8 #0.1
            conc=1.0
        else
            conc=0.0
        end
        return conc;

        #return 1.0;
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
            lon0=2.0*atan(1.0) #1.25*atan(1.0) #1.5*atan(1.0)
            d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
            if abs(d)<=0.1 #0.8 #0.1
                uS=1.0
            else
                uS=0.0
            end
            return velCa([uS,0.0,0.0],lon,lat)
        #end
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

    ha=zeros(Float64,m.topology.size[3])
    if case==:zPanel
        ha[[46,49,52]].=3.0
        ha[[47,50,53]].=2.0
        ha[[48,51,54]].=1.0
    elseif case==:xPanel
        ha[[19,22,25]].=3.0
        ha[[20,23,26]].=2.0
        ha[[21,24,27]].=1.0
    end
    p.solution[0.0].theta=ha

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"0", printSpherical=true)

    function V(xyz::Array{Float64,1})
      x=xyz[1]; y=xyz[2]; z=xyz[3];
      lon,lat,r=cart2sphere(x,y,z);
      uS=UMax*cos(lat)
      return velCa([uS,0.0,0.0],lon,lat)
    end
    Vfcomp=:RT0
    Vf=projectAdvection(p,V,Vfcomp);

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

    #RKadvection!(p, gamma, Vfcomp, Vf, nquadPhi, nquadPoints, MrT, MrV, 0.0, dt, EndTime, filename)

    y=p.solution[0.0];
    Y=Array{solution,1}(undef,MISMethod.nStage+1);
    FY=Array{solution,1}(undef,MISMethod.nStage);
    SthY=Array{SparseMatrixCSC{Float64,Int64},1}(undef,MISMethod.nStage);
    Time=0.0;

    for i=1:nIter
      @time y=splitExplicit(y,Y,FY,SthY,p,Vfcomp,Vf,gamma,nquadPhi,nquadPoints,MrT,MrV,MISMethod,Time,dt,ns);
      Time+=dt
      p.solution[Time]=y;
      p.solution[Time].theta=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoTheta,:rho,:rhoTheta,MrT);
      p.solution[Time].v=projectRhoChi(p,p.solution[Time].rho,p.solution[Time].rhoV,:rho,:rhoV,MrV)
      if mod(i,4)==0
          p2=deepcopy(p);
          unstructured_vtk(p2, Time, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"$i", printSpherical=true)
      end
      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename, printSpherical=true)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(0.0:200.0:EndTime)), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)

    return p;
end
