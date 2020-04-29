using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

gr()

include("solution.jl");
@solution(v,p,b);
include("constants.jl")
include("getKub.jl");
include("getQuad.jl");
include("meshTypes.jl");
include("meshFunctions.jl");
include("findall.jl");
include("transformation.jl")
include("initAssembledPhi.jl");

include("getOrderBoundary.jl")
include("getElementProperties.jl")
include("getCompoundElementProperties.jl")
include("degF.jl");
include("generateMesh.jl")
include("refineMesh.jl")
include("compoundData.jl")
include("getSubCells.jl")
include("assembleCompoundPhi.jl")
include("femProblem.jl");

include("l2g.jl")
include("jacobiCompound.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")
include("getDim.jl")
include("getSpace.jl")

include("assembMassCompound.jl");
include("assembLoadCompound.jl");
include("assembStiffCompound.jl");
include("applyStartValuesCompound.jl");

include("projectRecovery.jl")
include("embed.jl")
include("recovery.jl")
include("projectAdvection.jl")
include("advectionStiffCompound.jl")
include("discGalerkinCellsCompound.jl")
include("discGalerkinEdgesCompound.jl")
include("MIS.jl")
include("coordTrans.jl")
include("getEdgeType.jl")
include("setEdgeData.jl")
include("advectionBA.jl")
include("splitExplicitBA.jl")
include("symplektischerEulerBA.jl")

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtkCompound.jl");
include("vtkTest.jl");
include("jld.jl");
