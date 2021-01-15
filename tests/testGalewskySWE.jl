include("modulesSWE.jl")
#include("advectionStiffN.jl")

const stencilOrder=2;

const recoverySpace=Symbol("DGQuad")
const recoverySpaceVec=Symbol("VecDGQuadS")

@recovery(recoverySpace,recoverySpaceVec)

function testGalewsky()                                                                                                                     ^

    filename = "galewskyITR";
    nNodes=92
    #order: comp, compTest, recoverySpace

    femType=Dict(:h=>[:DG0],
                 :hV=>[:RT0, :RT0, recoverySpaceVec],
                 #:hV=>[:RT0, :VecP1S, :VecDG1S, :RT0B],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :f=>[:P1],);

    #=
    femType=Dict(:h=>[:DG0, :P1, :DG1, :DG0],
                 :hV=>[:RT0, :VecP1, :VecDG1, :RT0B],
                 :p=>[:DG0],
                 :v=>[:RT0],
                 :f=>[:P1],);

    #higher spaces

    femType=Dict(:h=>[:DG1, :P1, :DG1, :DG0],
                 :hV=>[:RT1, :VecP1, :VecDG1, :RT0B],
                 :p=>[:DG1],
                 :v=>[:RT1]);
    =#

    taskRecovery=true;
    adv=true;

    m=generateCubedSphere(92,6300000.0,0,:cube1)

    p=femProblem(m, femType,t=:shallow, advection=adv, taskRecovery=taskRecovery,
    stencilOrder=stencilOrder);

    gamma=0.5; #upwind
    MISMethod=MIS(:MIS2); #method of time integration

    dt=200.0  #160.0 #200.0 #50.0 #1200.0;
    ns=6;
    EndTime=6.0*24*60*60;
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
    function fh(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        return (Grav*H0G-(simpson(-0.5*pi,lat,r,pi/100.0,integrandG)))/Grav#+hH*cos(lat)*exp(-((lon-pi)/alphaG)^2.0)*exp(-((pi/4.0-lat)/betaG)^2.0)
    end
    uM=80.0
    lat0G=pi/7.0
    lat1G=pi/2.0-lat0G
    eN=exp(-4.0/(lat1G-lat0G)^2.0)
    function fvel(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        if lat<=lat0G || lat>=lat1G
          uS=0
        else
          uS=uM/eN*exp(1.0/((lat-lat0G)*(lat-lat1G)))
        end
        return velCa([uS,0.0,0.0],lon,lat)
    end
    f=Dict(:h=>fh,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

    h0=p.solution[0.0].h;
    p.solution[0.0].hV=projectChi(p,h0,p.solution[0.0].v,:h,:v);

    unstructured_vtk(p, 0.0, [:h, :hV, :v], ["h", "hV", "Velocity"], "testSphere/"*filename*"0", printSpherical=true)

    advectionTypes=Symbol[];
    for i in [:h,:hV]
        push!(advectionTypes,femType[i][1]);
        taskRecovery && push!(advectionTypes,femType[i][3]);
    end
    nquadPhi, nquadPoints=coordTrans(m, m.normals, advectionTypes, size(p.kubWeights,2));
    setEdgeData!(p, :v)
    recoveryMatrix!(p)

    MrV=assembMass(p.degFBoundary[femType[:hV][1]], m, p.kubPoints, p.kubWeights);

    y=p.solution[0.0];
    Y=Array{solution,1}(undef,MISMethod.nStage+1);
    FY=Array{solution,1}(undef,MISMethod.nStage);
    Time=0.0;

    nVTK=24*60*60/dt;
    #nVTK=2
    for i=1:nIter
        @time y=splitExplicit(y,Y,FY,p,gamma,nquadPhi,nquadPoints,MrV,MISMethod,Time,dt,ns);
        Time+=dt
        p.solution[Time]=y;
        p.solution[Time].v=projectRhoChi(p,p.solution[Time].h,p.solution[Time].hV,:h,:hV,MrV)
        if mod(i,nVTK)==0
            p2=deepcopy(p);
            unstructured_vtk(p2, Time, [:h, :hV, :v], ["h", "hV", "Velocity"], "testSphere/"*filename*"$i", printSpherical=false)
        end
        println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:h, :hV, :v], ["h", "hV", "Velocity"], "testSphere/"*filename, printSpherical=true)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)

    return p;
end
p=testGalewsky()
