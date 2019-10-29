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
include("solveB.jl");

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk.jl");
