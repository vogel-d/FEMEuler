using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

#pyplot()
gr() #Pkg.build("CodecZlib") #Pkg.build("LightXML")

include("solution.jl");
@solution(v,p,b);
include("constants.jl")
include("getKub.jl");
include("getQuad.jl");
include("meshTypes.jl");
include("meshFunctions.jl");
include("findall.jl");
include("transformation.jl")

include("coordTrans.jl")
include("compoundData.jl")

include("getOrderBoundary.jl")
include("getElementProperties.jl")
include("degF.jl");
include("generateMesh.jl")
include("refineMesh.jl")
include("femProblem.jl");
include("splitCompoundMesh.jl")
include("initAssembledPhi.jl")

include("l2g.jl")
include("jacobi.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")

include("assembMass.jl");
include("assembLoad.jl");
include("assembStiff.jl");
include("applyStartValues.jl");
include("solveB.jl");

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk.jl");
include("vtkTest.jl");
include("vtkRefined.jl");
include("jld.jl");
