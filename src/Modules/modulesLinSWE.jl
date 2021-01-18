using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;
using OffsetArrays;

gr()

include("../System/Variables/solution.jl");
@solution(h, v);
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
include("../FiniteElements/getSphereElementProperties.jl")
include("../FiniteElements/getElementProperties.jl")
include("../FiniteElements/degF.jl");
include("../Mesh/generateMesh.jl")
include("../Mesh/refineMesh.jl")
include("../System/femProblem.jl");

include("../FiniteElements/l2g.jl")
include("../Transformation/jacobiSphere.jl");
include("../Transformation/initJacobi.jl")
include("../Mesh/adaptGeometry.jl")
include("../Utilities/additionalFunctions.jl")

include("../System/Matrices/assembMass.jl");
include("../System/Matrices/assembLoad.jl");
include("../System/Matrices/assembStiff.jl");
include("../System/applyStartValues.jl");

include("../Output/plotSolution.jl");
include("../Output/plotSolutionGif.jl");
include("../Output/plotMesh.jl");
include("../Output/plotFEM.jl")
include("../Output/vtk.jl");
include("../Output/vtkTest.jl");
include("../Output/vtkSphere.jl");
include("../Output/vtkRefined.jl");
include("../Output/jld.jl");

include("../Mesh/Sphere/generateCubedSphere.jl")
include("../Mesh/Sphere/insert.jl")
include("../Mesh/Sphere/xmtoxc.jl")
include("../Mesh/Sphere/set_cuco.jl")
include("../Mesh/Sphere/tay.jl")
include("../Mesh/Sphere/fft.jl")
include("../Mesh/setOrientation.jl")
include("../Mesh/Sphere/cart2sphere.jl")
include("../Utilities/simpson.jl")
include("../Mesh/Sphere/vel.jl")
