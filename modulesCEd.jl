using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

gr()

include("solution.jl");
@solution(rho, rhoV, rhoTheta, v, theta);
include("diagnostic.jl");
@diagnostic(rhoBar, pBar, thBar);
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
include("femProblemd.jl");

include("l2g.jl")
include("jacobi.jl");
include("initJacobi.jl")
include("adaptGeometry.jl")
include("additionalFunctions.jl")

include("assembMass.jl");
include("assembMassRho.jl");
include("assembLoad.jl");
include("assembStiff.jl");
include("applyStartValuesd.jl");

include("projectRecovery.jl")
include("embed.jl")
include("recovery.jl")
include("advectionStiff.jl")
include("advectionStiffMatrix.jl")
include("discGalerkinCells.jl")
include("discGalerkinEdges.jl")
include("MIS.jl")
include("coordTrans.jl")
include("getEdgeType.jl")
include("setEdgeData.jl")
include("advectionCE.jl")
include("splitExplicitCE.jl")
include("symplektischerEulerCEd.jl")
include("projectPressure.jl")
include("projectChi.jl")
include("projectRhoChi.jl")


include("plotSolution.jl");
include("plotSolutionGif.jl");
include("plotMesh.jl");
include("plotFEM.jl")
include("vtk.jl");
include("jld.jl");
