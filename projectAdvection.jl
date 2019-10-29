function projectAdvection(p::femProblem,V::Array{Function,1},comp::Symbol)
    l=assembLoad(p.degFBoundary[comp],V,p.mesh,p.kubPoints,p.kubWeights);
    Vf=p.massMBoundary[comp]\l;
    return sparse(Vf);
end
