include("../src/Modules/modulesLinSWE.jl")

function testLinearSWE()

    filename = "linearSWE";

    #order: comp, compTest, recoverySpace
    femType=Dict(:h=>[:DG0], :v=>[:RT0], :f=>[:P1]);
    #higher spaces
    #femType=Dict(:h=>[:DG1],:v=>[:RT1]);

    taskRecovery=false;
    adv=false;

    n=30;
    m=generateCubedSphere(n,1.0,0,:cube1)

    p=femProblem(m, femType,t=:linshallow, advection=adv, taskRecovery=taskRecovery);

    dt=0.001;
    EndTime=10.0;
    nIter=Int64(EndTime/dt);

    #start functions
    function fh(xyz::Array{Float64,1})
        return exp(-((-xyz[1]-1.0)/0.1)^2);
    end
    function fvel(xyz::Array{Float64,1})
        return [0.0,0.0,0.0]
    end
    f=Dict(:h=>fh,:v=>fvel);

    assembMass!(p);
    assembStiff!(p);
    applyStartValues!(p, f);

    unstructured_vtk(p, 0.0, [:h, :v], ["h", "Velocity"], "testSphere/"*filename*"0", printSpherical=true)

    intMethod=:RungeKutta

    y=p.solution[0.0];
    Time=0.0;
    if intMethod==:Euler
        for i=1:nIter
          rh=p.massM[femType[:h][1]]\(p.stiffM[:div]*y.v);
          rVel=p.massM[femType[:v][1]]\(p.stiffM[:grad]*y.h);
          y.h+=dt*rh;
          y.v+=dt*rVel;
          Time+=dt
          p.solution[Time]=y;
          if mod(i,10)==0
              p2=deepcopy(p);
              unstructured_vtk(p2, Time, [:h, :v], ["h", "Velocity"], "testSphere/"*filename*"$i", printSpherical=false)
          end
          println(Time)
        end
    elseif intMethod==:RungeKutta
        for i=1:nIter
          rh=p.massM[femType[:h][1]]\(p.stiffM[:div]*y.v);
          rVel=p.massM[femType[:v][1]]\(p.stiffM[:grad]*y.h+p.stiffM[:coriolis]*y.v);
          hNeu=y.h+0.5*dt*rh;
          VelNeu=y.v+0.5*dt*rVel;

          rh=p.massM[femType[:h][1]]\(p.stiffM[:div]*VelNeu);
          rVel=p.massM[femType[:v][1]]\(p.stiffM[:grad]*hNeu+p.stiffM[:coriolis]*VelNeu);

          y.h+=dt*rh;
          y.v+=dt*rVel;

          Time+=dt
          p.solution[Time]=y;
          if mod(i,50)==0
              p2=deepcopy(p);
              unstructured_vtk(p2, Time, [:h, :v], ["h", "Velocity"], "testSphere/"*filename*"$i", printSpherical=false)
          end

          println(Time)
        end
    end

    #Speichern des Endzeitpunktes als vtu-Datei:
    unstructured_vtk(p, EndTime, [:h, :v], ["h", "Velocity"], "testSphere/"*filename, printSpherical=true)
    #Speichern aller berechneten Zwischenwerte als vtz-Datei:
    #unstructured_vtk(p, sort(collect(keys(p.solution))), [:h, :hV, :v], ["h", "hV", "Velocity"], "testSphere/"*filename)

    return p;
end
p=testLinearSWE();
