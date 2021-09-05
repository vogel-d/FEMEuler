include("../src/Modules/modulesCEAdv.jl")

function testAdvS()
    filename = "AdvS";

    case=:periodic
    stencilOrder=1;
    recoveryOrder=1;

    recoverySpace=Symbol("R$recoveryOrder")

    #order: comp, compHigh, compRec, compDG
    femType=Dict(:rho=>[:DG0, :DG0, recoverySpace],
                 :rhoV=>[:RT0, :RT0, :VecDG1],
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

    m=generateRectMesh(4,4,:periodic,:periodic); #(east/west, top/bottom)

    #adaptGeometry!(m,(0.3,0.3),false); #sin perbutation

    p=femProblem(m, femType, t=:compressible, advection=advection,
        taskRecovery=taskRecovery, stencilOrder=stencilOrder, g=2);

    gamma=0.5; #upwind
    UMax=20.0; #UMax determines the advection in x direction
    MISMethod=MIS(:MIS2); #method of time integration

    dt=2.0;
    ns=15;
    EndTime=1000.0;
    nIter=Int64(EndTime/dt);

    #start functions
    function frho(xz::Array{Float64,1})
        return 1.0;
    end
    function ftheta(xz::Array{Float64,1})
        return 0.0;
    end
    function fvel(xz::Array{Float64,1})
        return [UMax, 0.0]
    end
    f=Dict(:rho=>frho,:theta=>ftheta,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    #p.boundaryValues[(:theta,:P1)]=300.0*ones(p.degFBoundary[:P1].numB-p.degFBoundary[:P1].num);
    applyStartValues!(p, f);

    ha=zeros(Float64,m.topology.size[3])
    nx=m.topology.n[1]; ny=m.topology.n[2];
    if case==:linx
        #val=nx/2.0
        val=Float64(nx)
        for i in 1:nx
            ha[i:nx:nx*(ny-1)+i].=val
            #if iseven(i)
                val-=1.0
            #end
        end

    elseif case==:linz
        val=ny/2.0-1.0
        z=1;
        for i in 1:ny
            ha[z:(nx-1)+z].=val
            z+=nx;
            if iseven(i)
                val-=1.0
            end
        end
    elseif case==:quadx
        val=Float64(nx)
        for i in 1:nx
            ha[i:nx:nx*(ny-1)+i].=val^2
            val-=1.0
        end
    elseif case==:linx2
        xmin=-5.0
        xmax=5.0
        val=Float64(nx)
        for i in 1:nx
            ha[i:nx:nx*(ny-1)+i].=(val-xmin)/(xmax-xmin)
            val-=1.0
        end
    elseif case==:periodic
        #Für 4x4-Gitter

        ha[[1,5,9,13]].=2.0
        ha[[2,6,10,14]].=1.0
        ha[[3,7,11,15]].=4.0
        ha[[4,8,12,16]].=3.0

        #=
        ha[[1,5,9,13]].=3.0
        ha[[2,6,10,14]].=2.0
        ha[[3,7,11,15]].=1.0
        ha[[4,8,12,16]].=4.0
        =#
    elseif case==:max
        #Für 8x8-Gitter
        ha[[28,29,36,37]].=1.0
    end
    p.solution[0.0].theta=ha

    rho0=p.solution[0.0].rho;
    p.solution[0.0].rhoTheta=projectChi(p,rho0,p.solution[0.0].theta,:rho,:theta);
    p.solution[0.0].rhoV=projectChi(p,rho0,p.solution[0.0].v,:rho,:v);

    unstructured_vtk(p, 0.0, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"0")

    function V(xz::Array{Float64,1})
      return [UMax, 0.0]
    end
    Vfcomp=:RT0
    Vf=projectAdvection(p,V,Vfcomp);

    #taskRecovery ? pos=[1,3] : pos=[1];
    advectionTypes=Symbol[];
    recoveryTypes=Symbol[];
    for i in [:rho,:rhoTheta,:rhoV]
        push!(advectionTypes,femType[i][1]);
        if i==:rhoV
            taskRecovery && push!(advectionTypes,femType[i][3]);
        else
            taskRecovery && push!(recoveryTypes,femType[i][3]);
        end
    end
    nquadPhi, nquadPoints=coordTrans(m, m.normals, advectionTypes, recoveryTypes, size(p.kubWeights,2));
    setEdgeData!(p, :v)

    MrT=assembMass(p.degFBoundary[femType[:rhoTheta][1]], m, p.kubPoints, p.kubWeights);
    MrV=assembMass(p.degFBoundary[femType[:rhoV][1]], m, p.kubPoints, p.kubWeights);

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
          unstructured_vtk(p2, Time, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename*"$i")
      end
      println(Time)
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:rho, :rhoV, :rhoTheta, :v, :theta], ["Rho", "RhoV", "RhoTheta", "Velocity", "Theta"], "testCompressibleEuler/"*filename)

    return p
end
