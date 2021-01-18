using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

gr()

include("../System/Variables/solution.jl");
@solution(rho, rhoV, rhoTheta, v, theta);
include("../System/constants.jl")
include("../Utilities/getKub.jl");
include("../Utilities/getQuad.jl");
include("../Mesh/meshTypes.jl");
include("../Mesh/meshFunctions.jl");
include("../Utilities/findall.jl");
include("../Transformation/transformation.jl")

include("../System/Compound/compoundData.jl")
include("../System/femData.jl")
include("../System/getOrderBoundary.jl")
include("../FiniteElements/getQuadElementProperties.jl")
include("../FiniteElements/getSphereElementProperties.jl")
include("../FiniteElements/getTriElementProperties.jl")
include("../FiniteElements/getElementProperties.jl")
include("../FiniteElements/degF.jl");
include("../Mesh/generateMesh.jl")
include("../Mesh/refineMesh.jl")
include("../System/femProblem.jl");
include("../Mesh/splitCompoundMesh.jl")
include("../System/Compound/initAssembledPhi.jl")

include("../FiniteElements/l2g.jl")
include("../Transformation/jacobi.jl");
include("../Transformation/initJacobi.jl")
include("../Mesh/adaptGeometry.jl")
include("../Utilities/additionalFunctions.jl")

include("../System/Matrices/assembMass.jl");
include("../System/Matrices/assembMassRho.jl");
include("../System/Matrices/assembLoad.jl");
include("../System/Matrices/assembStiff.jl");
include("../System/applyStartValues.jl");

include("../Advection/Recovery/projectRecovery.jl")
include("../Advection/Recovery/embed.jl")
include("../Advection/Recovery/recovery.jl")
include("../Advection/AdvectionStiff/advectionStiff.jl")
include("../Advection/AdvectionStiff/advectionStiffMatrix.jl")
include("../Advection/DiscontinousGalerkin/discGalerkinCells.jl")
include("../Advection/DiscontinousGalerkin/discGalerkinEdges.jl")
include("../Solvers/MIS.jl")
include("../Transformation/coordTrans.jl")
include("../Advection/setEdgeData.jl")
include("../Advection/advectionCE.jl")
include("../Solvers/splitExplicitCE.jl")
include("../Solvers/symplektischerEulerCE.jl")
include("../System/Variables/projectPressure.jl")
include("../System/Variables/projectChi.jl")
include("../System/Variables/projectRhoChi.jl")

include("../Output/plotSolution.jl");
include("../Output/plotSolutionGif.jl");
include("../Output/plotMesh.jl");
include("../Output/plotFEM.jl")
include("../Output/vtk.jl");
include("../Output/vtkTest.jl");
include("../Output/vtkRefined.jl");
include("../Output/jld.jl");

include("../Advection/Recovery/getStencil.jl")
include("../Advection/Recovery/assembRecovery.jl")
include("../FiniteElements/getPhiRecovery.jl")

include("../Advection/AdvectionStiff/advectionStiffR.jl")
include("../Advection/DiscontinousGalerkin/discGalerkinCellsR.jl")
include("../Advection/DiscontinousGalerkin/discGalerkinEdgesR.jl")
include("../Output/vtkRecovery.jl");
include("../Advection/Recovery/transformRecoveryCoord.jl")
include("../Advection/Recovery/getTangentialPlane.jl")
include("../Advection/Recovery/recoveryMatrix.jl")
include("../Advection/Recovery/intersectPlane.jl")
