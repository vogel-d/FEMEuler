using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;
#using Distributed;
using OffsetArrays;

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
include("jacobiSphere.jl");
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
include("advectionCE.jl")
include("splitExplicitCE.jl")
include("symplektischerEulerCE.jl")
include("projectPressure.jl")
include("projectChi.jl")
include("projectRhoChi.jl")


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

include("getStencil.jl")
include("assembRecovery.jl")
include("getPhiRecovery.jl")

include("advectionStiffRecovery.jl")
#include("advectionStiffMatrix.jl")
include("discGalerkinCellsRecovery.jl")
include("discGalerkinEdgesRecovery.jl")
include("vtkRecovery.jl");
include("transformRecoveryCoord.jl")
include("getTangentialPlane.jl")
include("recoveryMatrix.jl")
