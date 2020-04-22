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
include("getSubCells.jl")
include("getSpace.jl")
include("transformation.jl")

include("getOrderBoundary.jl")
include("getElementProperties.jl")
include("degF.jl");
include("generateMesh.jl")
include("refineMesh.jl")
include("assembleCompoundPhi.jl")
include("compoundData.jl")
include("femProblem.jl");

include("l2g.jl")
include("jacobiCompound.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")

include("assembMassCompound.jl");
include("assembLoadCompound.jl");
include("assembStiffCompound.jl");
include("applyStartValuesCompound.jl");
include("solveB.jl");

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk.jl");
include("vtkTest.jl");
include("jld.jl");
