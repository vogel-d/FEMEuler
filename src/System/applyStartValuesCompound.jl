function applyStartValuesCompound!(p::femProblem,f)
    sol=solution();
    for i in fieldnames(solution)
        nc=p.degFBoundary[p.femType[i][1]].numB
        ni=p.degFBoundary[p.femType[i][1]].num
        nb=nc-ni;
        ti=p.femType[i][1];
        if haskey(f,i)
            F=p.massMBoundary[ti];
            l=assembLoadCompound(p.degFBoundary[ti],f[i], p.mesh, p.kubPoints, p.kubWeights, p.data.compoundData)
            h=F\l;
        else
            h=zeros(nc)
        end
        if iszero(nb)
            setfield!(sol,i,h);
        elseif haskey(p.boundaryValues,(i,ti))
            h[(ni+1):nc]=p.boundaryValues[(i,ti)]
            setfield!(sol,i,h);
        else
            h[(ni+1):nc]=zeros(nb);
            setfield!(sol,i,h);
        end
    end
    p.solution[0.0]=sol;
    return nothing;
end
