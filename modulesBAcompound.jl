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
include("degF.jl");
include("generateMesh.jl")
include("refineMesh.jl")
include("compoundData.jl")
include("getCompoundElementProperties.jl")
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
include("orderAdjacentSubCells.jl");

include("assembMassCompound.jl");
include("assembLoadCompound.jl");
include("assembStiffCompound.jl");
include("applyStartValuesCompound.jl");

include("projectRecovery.jl")
include("embed.jl")
include("recovery.jl")
include("projectAdvectionCompound.jl")
include("advectionStiffCompound.jl")
include("advectionStiffCompoundMatrix.jl")
include("discGalerkinCellsCompound.jl")
include("discGalerkinEdgesCompound.jl")
include("MIS.jl")
include("coordTrans.jl")
include("getEdgeType.jl")
include("correctNormalsCompound.jl")
include("setCompoundEdgeData.jl")
include("advectionBA.jl")
include("splitExplicitBA.jl")
include("symplektischerEulerBA.jl")

include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("splitCompoundMesh.jl")
include("vtkCompound.jl");
include("vtkTest.jl");
include("jld.jl");