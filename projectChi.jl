#Aus chi und rho wird rhoChi bestimmt
function projectChi(p::femProblem,valRho::Array{AbstractFloat,1},valChi::Array{AbstractFloat,1},compRho::Symbol, compChi::Symbol)
    compChi=p.femType[compChi][1];
    MRho=assembMassRho(p.degFBoundary[compChi], p.degFBoundary[p.femType[compRho][1]], valRho, p.mesh, p.kubPoints, p.kubWeights);
    return p.massMBoundary[compChi]\(MRho*valChi);
end
