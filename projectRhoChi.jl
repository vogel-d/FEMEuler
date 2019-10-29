#Aus rhoChi und rho wird chi bestimmt
function projectRhoChi(p::femProblem,valRho::Array{Float64,1},valChi::Array{Float64,1},compRho::Symbol, compChi::Symbol,M::SparseMatrixCSC{Float64,Int64})
    compChi=p.femType[compChi][1];
    FRho=lu(assembMassRho(p.degFBoundary[compChi], p.degFBoundary[p.femType[compRho][1]], valRho, p.mesh, p.kubPoints, p.kubWeights));
    #M=assembMass(p.degFBoundary[compChi], p.mesh, p.kubPoints, p.kubWeights);
    return FRho\(M*valChi);
end
