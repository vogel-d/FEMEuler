using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

gr()

include("solution.jl");
@solution(rho, rhoV, rhoTheta, v, theta);
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
include("jacobi.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")

include("assembMass.jl");
include("assembMassRho.jl");
include("assembLoad.jl");
include("assembStiff.jl");
include("applyStartValues.jl");

include("projectRecovery.jl")
include("embed.jl")
include("recovery.jl")
include("advectionStiff.jl")
include("advectionStiffMatrix.jl")
include("discGalerkinCells.jl")
include("discGalerkinEdges.jl")
include("MIS.jl")
include("coordTrans.jl")
include("setEdgeData.jl")
include("advectionCEAdv.jl")
include("splitExplicitCEAdv.jl")
include("symplektischerEulerCEAdv.jl")
include("projectPressure.jl")
include("projectChi.jl")
include("projectRhoChi.jl")


include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk.jl");
include("vtkTest.jl");
include("vtkRefined.jl");
include("jld.jl");

include("projectAdvection.jl")

include("getStencil.jl")
include("assembRecovery.jl")
include("getPhiRecovery.jl")

include("advectionStiffR.jl")
include("discGalerkinCellsR.jl")
include("discGalerkinEdgesR.jl")
include("vtkRecovery.jl");
include("transformRecoveryCoord.jl")
include("getTangentialPlane.jl")
include("intersectPlane.jl")
