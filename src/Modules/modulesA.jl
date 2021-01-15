using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

gr()

include("../System/Variables/solution.jl");
@solution(v,p,b);
include("../Utilities/getKub.jl");
include("../Utilities/getQuad.jl");
include("../Mesh/meshTypes.jl");
include("../Mesh/meshFunctions.jl");
include("../Utilities/findall.jl");
include("../Transformation/transformation.jl")

include("../System/getOrderBoundary.jl")
include("../FiniteElements/getQuadElementProperties.jl")
include("../FiniteElements/getSphereElementProperties.jl")
include("../FiniteElements/getTriElementProperties.jl")
include("../FiniteElements/getElementProperties.jl")
include("../FiniteElements/degF.jl");
include("../Mesh/generateMesh.jl")
include("../Mesh/refineMesh.jl")
include("../System/Compound/compoundData.jl")
include("../System/femData.jl")
include("../System/femProblem.jl");

include("../FiniteElements/l2g.jl")
include("../Transformation/jacobi.jl");
include("../Transformation/initJacobi.jl")
include("../Mesh/adaptGeometry.jl")
include("../Utilities/additionalFunctions.jl")

include("../System/Matrices/assembMass.jl");
include("../System/Matrices/assembLoad.jl");
include("../System/Matrices/assembStiff.jl");
include("../System/applyStartValues.jl");

include("../Advection/Recovery/projectRecovery.jl")
include("../Advection/Recovery/embed.jl")
include("../Advection/Recovery/recovery.jl")
include("../Advection/projectAdvection.jl")
include("../Advection/AdvectionStiff/advectionStiff.jl")
include("../Advection/DiscontinousGalerkin/discGalerkinCells.jl")
include("../Advection/DiscontinousGalerkin/discGalerkinEdges.jl")
include("../Solvers/MIS.jl")
include("../Transformation/coordTrans.jl")
include("../Advection/setEdgeData.jl")
include("../Advection/advectionBA.jl")
include("../Solvers/splitExplicitBA.jl")
include("../Solvers/symplektischerEulerA.jl")
include("../Solvers/RKadvectionCE.jl")

include("../Output/plotSolution.jl");
include("../Output/plotSolutionGif.jl");
include("../Output/plotMesh.jl");
include("../Output/plotFEM.jl")
include("../Output/vtk.jl");
include("../Output/vtkTest.jl");
include("../Output/jld.jl");
