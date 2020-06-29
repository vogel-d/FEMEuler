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
include("getSpace.jl")
include("transformation.jl")
include("initAssembledPhi.jl");
include("coordTrans.jl");

include("getOrderBoundary.jl")
include("getElementProperties.jl")
include("degF.jl");
include("generateMesh.jl")
include("refineMesh.jl")
include("compoundData.jl")
include("getCompoundElementProperties.jl")
include("getSubCells.jl")
include("femProblem.jl");
include("assembleCompoundPhi.jl")
include("assemblePhiPre.jl")

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
include("splitCompoundMesh.jl")
include("vtkCompound.jl");
include("vtkTest.jl");
include("jld.jl");
