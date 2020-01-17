function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1}) #,m::mesh,kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    m=p.mesh;
    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    degF=p.degFBoundary;
    cH=projectRecovery(degF[comp[2]],degF[comp[1]],cval,p.massMBoundary[comp[2]],m,kubPoints,kubWeights);
    cHP=projectRecovery(degF[comp[4]],degF[comp[2]],cH,p.massMBoundary[comp[4]],m,kubPoints,kubWeights);

    n=m.topology.size[3]

    cR=embed(comp[2],degF[comp[2]],cH,comp[3],degF[comp[3]],n);
    cEmbed=embed(comp[1],degF[comp[1]],cval,comp[3],degF[comp[3]],n);
    cHPEmbed=embed(comp[4],degF[comp[4]],cHP,comp[3],degF[comp[3]],n);
    
    return cR+(cEmbed-cHPEmbed);
end

function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1}, bcomp::Symbol)
    m=p.mesh;
    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    degF=p.degFBoundary;

    cH=projectRecovery(degF[comp[2]],degF[comp[1]],cval,p.massMBoundary[comp[2]],m,kubPoints,kubWeights);
    if haskey(p.boundaryValues,(bcomp,comp[2]))
        cH[degF[comp[2]].num+1:end]=p.boundaryValues[(bcomp,comp[2])];
    end
    cHP=projectRecovery(degF[comp[4]],degF[comp[2]],cH,p.massMBoundary[comp[4]],m,kubPoints,kubWeights);
    if haskey(p.boundaryValues,(bcomp,comp[4]))
        cH[degF[comp[4]].num+1:end]=p.boundaryValues[(bcomp,comp[4])];
    end

    n=m.topology.size[3]

    cR=embed(comp[2],degF[comp[2]],cH,comp[3],degF[comp[3]],n);
    cEmbed=embed(comp[1],degF[comp[1]],cval,comp[3],degF[comp[3]],n);
    cHPEmbed=embed(comp[4],degF[comp[4]],cHP,comp[3],degF[comp[3]],n);

    return cR+(cEmbed-cHPEmbed);
end
