using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;
using Distributed;

pyplot()

include("solution.jl");
@solution(rho, rhoV, rhoTheta, v, theta);
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
include("correctVelocity.jl")


include("generateEquals.jl")
include("generateMixedBoundary.jl")
include("generatePeriodicBoundary.jl")
include("assembMass.jl");
include("assembMassRho.jl");
include("assembLoad.jl");
include("assembStiff.jl");
include("assembFEM.jl");
include("applyStartValues.jl");

include("projectRecovery.jl")
include("embed.jl")
include("recovery.jl")
include("advectionStiff.jl")
include("advectionStiffMatrix.jl")
include("MIS.jl")
include("coordTrans.jl")
include("getPhi.jl")
include("getEdgeType.jl")
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
include("vtk.jl");
