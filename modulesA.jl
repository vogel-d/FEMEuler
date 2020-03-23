using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

gr()

include("solution.jl");
@solution(v,p,b);
include("getKub.jl");
include("getQuad.jl");
include("meshTypes.jl");
include("meshFunctions.jl");
include("findall.jl");
include("transformation.jl")

include("getOrderBoundary.jl")
include("getElementProperties.jl")
include("degF.jl");
include("generateMesh.jl")
include("refineMesh.jl")
include("femProblem.jl");

include("l2g.jl")
include("jacobi.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")
include("getDim.jl")
include("getSpace.jl")

include("assembMass.jl");
include("assembLoad.jl");
include("assembStiff.jl");
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
include("getEdgeType.jl")
include("setEdgeData.jl")
include("advectionA.jl")
include("splitExplicitBA.jl")
include("symplektischerEulerA.jl")
include("RKadvection.jl")

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk2D.jl");
include("jld.jl");
