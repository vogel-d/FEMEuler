function projectAdvectionCompound(p::femProblem,V::Array{Function,1},comp::Symbol)
    l=assembLoadCompound(p.degFBoundary[comp],V,p.mesh,p.kubPoints,p.kubWeights, p.data.compoundData);
    Vf=p.massMBoundary[comp]\l;
    return sparse(Vf);
end
