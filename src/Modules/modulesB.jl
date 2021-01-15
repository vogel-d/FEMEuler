using SparseArrays;
using LinearAlgebra;
using SuiteSparse;
using Plots;

#pyplot()
gr() #Pkg.build("CodecZlib") #Pkg.build("LightXML")

include("../System/Variables/solution.jl");
@solution(v,p,b);
include("../System/constants.jl")
include("../Utilities/getKub.jl");
include("../Utilities/getQuad.jl");
include("../Mesh/meshTypes.jl");
include("../Mesh/meshFunctions.jl");
include("../Utilities/findall.jl");
include("../Transformation/transformation.jl")

include("../Transformation/coordTrans.jl")
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
include("../System/Matrices/assembLoad.jl");
include("../System/Matrices/assembStiff.jl");
include("../System/applyStartValues.jl");
include("../Solvers/solveB.jl");

include("../Output/plotSolution.jl");
include("../Output/plotSolutionGif.jl");
include("../Output/plotMesh.jl");
include("../Output/plotFEM.jl")
include("../Output/vtk.jl");
include("../Output/vtkTest.jl");
include("../Output/vtkRefined.jl");
include("../Output/jld.jl");
