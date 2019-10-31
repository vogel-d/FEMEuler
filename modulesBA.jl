using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;
using Distributed;

pyplot()

include("solution.jl");
@solution(v,p,b);
include("getKub.jl");
include("getQuad.jl");
include("meshTypes.jl");
include("meshFunctions.jl");
include("findall.jl");
include("transformation.jl")

include("getElementProperties.jl")
include("degF.jl");
include("getBoundary.jl")
include("generateMesh.jl")
include("femProblem.jl");

include("l2g.jl")
include("jacobi.jl");
include("initPhi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")

include("generateEquals.jl")
include("generateMixedBoundary.jl")
include("generatePeriodicBoundary.jl")
include("assembMass.jl");
include("assembLoad.jl");
include("assembStiff.jl");
include("assembFEM.jl");
include("applyStartValues.jl");

include("projectRecovery.jl")
include("embed.jl")
include("recovery.jl")
include("projectAdvection.jl")
include("advectionStiff.jl")
include("discGalerkinCells.jl")
include("discGalerkinEdges.jl")
include("MIS.jl")
include("coordTrans.jl")
include("getPhi.jl")
include("getEdgeType.jl")
include("setEdgeData.jl")
include("advectionBA.jl")
include("splitExplicitBA.jl")
include("symplektischerEulerBA.jl")

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk.jl");
