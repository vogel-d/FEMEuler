include("modulesSphereAdv.jl")
#include("advectionStiffN.jl")

stencilOrder=2;
recoveryOrder=1;

#recoverySpace=Symbol("R$recoveryOrder")
#recoverySpaceVec=Symbol("VecR$(recoveryOrder)S")

recoverySpace=Symbol("DGQuad")
recoverySpaceVec=Symbol("VecDG$(recoveryOrder)")

@recovery(recoverySpace,recoverySpaceVec)

function testSphereAdv2()

    filename = "testSphereVecAdvRN";

    #order: comp, compTest, recoverySpace
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 :rhoV=>[:RT0, :RT0, recoverySpaceVec],
                 #:rhoV=>[:RT0, :VecP1, :VecDG1, :RT0B],
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
    n=10
    Rad=1.0
    UMax=1.0
    #m=generateCubedSphere(40,6300000.0,0,:cube1)
    #m=generateCubedSphere(40,6300000.0)
    #m=generateCubedSphere(15,6300000.0,0,:cube1)
    m=generateCubedSphere(n,Rad,0,:cube1)
    #m=generateCubedSphere(3,6300000.0,0,:cube1)

    p=femProblem(m, femType,t=:compressible, advection=adv,
        taskRecovery=taskRecovery, recoveryOrder=recoveryOrder,stencilOrder=stencilOrder);

    gamma=0.5; #upwind
    #UMax=100.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MISRK2); #method of time integration

    CFL=0.25
    dt=CFL*2*pi/(4*n)/UMax*Rad
    ns=15;
    EndTime=80000.0
    #nIter=Int64(EndTime/dt);
    nIter=200

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
        if abs(d)<=0.8
        #if abs(d)<=0.4 #0.8 #0.1
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
            lon0=2.0*atan(1.0) #1.25*atan(1.0) #1.5*atan(1.0)
            d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0))
            if abs(d)<=0.8 #0.8 #0.1
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

    ha=zeros(Float64,m.topology.size[3])
    #ha[5]=1.0
    #ha[14]=1.0
    ha[23]=1.0
    #ha[32]=1.0
    #=
    for i in 1:9
        ha[i]=Float64(i)
    end
    =#
    #p.solution[0.0].theta=ha


    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    #refined_vtk(p, 2, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"RF", printSpherical=true)

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename*"0", printSpherical=true)

    function V(xyz::Array{Float64,1})
      x=xyz[1]; y=xyz[2]; z=xyz[3];
      return [-y,x,0.0]
    end

    Vfcomp=:RT0
    Vf=projectAdvection(p,V,Vfcomp);
    vtk(m,p.degFBoundary[Vfcomp],Vector(Vf),Vfcomp,"testVf", printSpherical=false)

    advectionTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        push!(advectionTypes,femType[i][1]);
        #(taskRecovery && length(femType[i])==4) && push!(advectionTypes,femType[i][3]);
        push!(advectionTypes,femType[i][3]);
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
    #=
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
    end
    =#

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
      println("Time ",Time," dt ",dt)
      if Time>=0.5*pi
          break
      end
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    #unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename, printSpherical=true)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(0.0:200.0:EndTime)), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["h", "hV", "hTheta", "Velocity", "Theta"], "testSphere/"*filename)

    return p;
end
