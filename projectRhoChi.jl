#Aus rhoChi und rho wird chi bestimmt
function projectRhoChi(p::femProblem,valRho::Array{AbstractFloat,1},valChi::Array{AbstractFloat,1},compRho::Symbol, compChi::Symbol,M::SparseMatrixCSC{AbstractFloat,Int})
    FRho=lu(assembMassRho(p.degFBoundary[p.femType[compChi][1]], p.degFBoundary[p.femType[compRho][1]], valRho, p.mesh, p.kubPoints, p.kubWeights));
    #M=assembMass(p.degFBoundary[compChi], p.mesh, p.kubPoints, p.kubWeights);
    return FRho\(M*valChi);
end
