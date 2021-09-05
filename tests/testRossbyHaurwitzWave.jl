include("../src/Modules/modulesSWE.jl")

const stencilOrder=2;

const recoverySpace=Symbol("DGQuad")
const recoverySpaceVec=Symbol("VecDGQuadS")

@recovery(recoverySpace,recoverySpaceVec)

function testRossbyHauritzWave()

    filename = "rossbyHauritzWaveRQuad50n";
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

    m=generateCubedSphere(50,6.37122e6,0,:cube1)

    p=femProblem(m, femType,t=:shallow, advection=adv, taskRecovery=taskRecovery,
    stencilOrder=stencilOrder);

    gamma=0.5; #upwind
    MISMethod=MIS(:MIS2); #method of time integration

    dt=20.0;
    ns=6;
    EndTime=14.0*24*60*60;
    nIter=Int64(EndTime/dt);

    #start functions
    H06=8000.0 #m
    omega6=7.8480e-6 #s^-1
    K6=7.8480e-6 #s^-1
    R6=4.0
    function fh(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        A=0.5*omega6*(2.0*Omega+omega6)*cos(lat)*cos(lat)+
            0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat)+
            (2.0*R6*R6-R6-2.0)-2.0*R6*R6*(cos(lat))^(-2))
        B=2.0*(Omega+omega6)*K6/(R6+1.0)/(R6+2.0)*
            (cos(lat))^R6*((R6*R6+2.0*R6+2.0)-
            (R6+1.0)*(R6+1.0)*cos(lat)*cos(lat))
        C=0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat)-(R6+2.0))
        return (Grav*H06+r*r*(A+B*cos(R6*lon)+C*cos(2.0*R6*lon)))/Grav
    end
    function fvel(xyz::Array{Float64,1})
        x=xyz[1]; y=xyz[2]; z=xyz[3];
        lon,lat,r=cart2sphere(x,y,z);
        uS=r*omega6*cos(lat)+
            r*K6*(cos(lat))^(R6-1.0)*
            (R6*sin(lat)*sin(lat)-cos(lat)*cos(lat))*cos(R6*lon)
        vS=-r*K6*R6*(cos(lat))^(R6-1.0)*sin(lat)*sin(R6*lon)
        velX=(-sin(lon)*cos(lat)*uS-cos(lon)*sin(lat)*vS)
        velY=(cos(lon)*cos(lat)*uS-sin(lon)*sin(lat)*vS)
        velZ=vS*cos(lat)
        #return velCa([velX,velY,velZ],lon,lat)
        return [velX,velY,velZ]
    end
    f=Dict(:h=>fh,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

    h0=p.solution[0.0].h;
    p.solution[0.0].hV=projectChi(p,h0,p.solution[0.0].v,:h,:v);

    unstructured_vtk(p, 0.0, [:h, :hV, :v], ["h", "hV", "Velocity"], "testSphere/"*filename*"0", printSpherical=true)

    advectionTypes=Symbol[];
    for i in [:hV]
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

    nVTK=6*60*60/dt;
    #nVTK=10
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
p=testRossbyHauritzWave();
