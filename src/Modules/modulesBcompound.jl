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
include("../System/Compound/initAssembledPhi.jl");
include("../Transformation/coordTrans.jl");

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
include("../FiniteElements/getCompoundElementProperties.jl")
include("../System/Compound/getSubCells.jl")
include("../System/femProblem.jl");
include("../System/Compound/assembCompoundPhi.jl")
include("../System/Compound/assembPhiPre.jl")

include("../FiniteElements/l2g.jl")
include("../Transformation/jacobiCompound.jl");
include("../Transformation/initJacobi.jl")
include("../Mesh/adaptGeometry.jl")
include("../Utilities/additionalFunctions.jl")

include("../System/Matrices/assembMassCompound.jl");
include("../System/Matrices/assembLoadCompound.jl");
include("../System/Matrices/assembStiffCompound.jl");
include("../System/applyStartValuesCompound.jl");
include("../Solvers/solveB.jl");

include("../Output/plotSolution.jl");
include("../Output/plotSolutionGif.jl");
include("../Output/plotMesh.jl");
include("../Output/plotFEM.jl")
include("../Mesh/splitCompoundMesh.jl")
include("../Output/vtkCompound.jl");
include("../Output/vtkTest.jl");
include("../Output/jld.jl");
