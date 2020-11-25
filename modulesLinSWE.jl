using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;
using OffsetArrays;

gr()

include("solution.jl");
@solution(h, v);
include("constants.jl")
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
include("jacobiSphere.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")

include("assembMass.jl");
include("assembLoad.jl");
include("assembStiff.jl");
include("applyStartValues.jl");

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("jld.jl");

include("generateCubedSphere.jl")
include("insert.jl")
include("xmtoxc.jl")
include("set_cuco.jl")
include("tay.jl")
include("fft.jl")
include("vtk.jl");
include("vtkTest.jl");
include("vtkSphere.jl");
include("vtkRefined.jl");
include("setOrientation.jl")
include("cart2sphere.jl")
include("simpson.jl")
include("vel.jl")
