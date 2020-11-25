function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1})
    return recovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[1]],p.degFBoundary[comp[3]],cval,p.stencil,p.recoveryM[(comp[2],comp[3])],p.mesh,p.kubPoints,p.kubWeights)
end

function recovery(degFT::degF{1,:H1}, degFF::degF{1,:H1}, degFR::degF{1,:H1}, cval::Array{Float64,1},
                  stencil::Array{Array{Int,1},1}, recoveryM::Array{Any,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})

    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    phiR=@views degFR.phi;

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    nT=length(phiT)
    nF=length(phiF)
    nR=length(phiR)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);
    globalNumR=Array{Int64,1}(undef,nR);
    cS=zeros(degFT.numB)
    cR=zeros(degFR.numB);
    for f in 1:nf
        cS[recoveryM[f][2]].=0.0
        for k in stencil[f]
            jacobi!(J,dJ,m,k,kubPoints,coord);
            k==f ? w=10.0^8 : w=1.0;

            l2g!(globalNumT,degFT,k)
            l2g!(globalNumF,degFF,k);

            for i in 1:nT
                gi=globalNumT[i];
                z=0.0;
                for j in 1:nF
                    gj=globalNumF[j];
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            z+=phiT[i][l,r]*kubWeights[l,r]*abs(dJ[l,r])*cval[gj]*phiF[j][l,r];
                        end
                    end
                end
                cS[gi] += w.*z;
            end
        end
        l2g!(globalNumR,degFR,f)
        cR[globalNumR]=recoveryM[f][1]\cS[recoveryM[f][2]]
    end
    return cR;
end

function recovery(degFT::degF{2,:H1div}, degFF::degF{2,:H1div}, degFR::degF{2,:H1xH1}, cval::Array{Float64,1},
                  stencil::Array{Array{Int,1},1}, recoveryM::Array{Any,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    phiR=@views degFR.phi;

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    nT=size(phiT,2)
    nF=size(phiF,2)
    nR=size(phiR,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    jphiF=initJacobi((m.geometry.dim,nF),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    rhs=Float64[]

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);
    globalNumR=Array{Int64,1}(undef,nR);
    cS=zeros(degFT.numB)
    cR=zeros(degFR.numB);
    for f in 1:nf
        cS[recoveryM[f][2]].=0.0
        for k in stencil[f]
            jacobi!(J,ddJ,jphiT,jphiF,m,k,kubPoints,phiT,phiF,coord);
            k==f ? w=10.0^8 : w=1.0;

            l2g!(globalNumT,degFT,k)
            l2g!(globalNumF,degFF,k);

            for i in 1:nT
                gi=globalNumT[i];
                z=0.0;
                for j in 1:nF
                    gj=globalNumF[j];
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            vecdot=0.0
                            for d in 1:m.geometry.dim
                                vecdot+=jphiT[d,i][l,r]*jphiF[d,j][l,r]
                            end
                            z+=kubWeights[l,r]*ddJ[l,r]^2/abs(ddJ[l,r])*cval[gj]*vecdot;
                        end
                    end
                end
                cS[gi] += w.*z;
            end
        end
        if m.geometry.dim==2
            rhs=cS[recoveryM[f][2]]
        else
            rhs=[cS[recoveryM[f][2]];zeros(length(stencil[f]))]
        end
        l2g!(globalNumR,degFR,f)
        cR[globalNumR]=recoveryM[f][1]\rhs
    end
    return cR;
end


function recovery(p::femProblem, dim::Int64, comp::Array{Symbol,1}, cval::Array{Float64,1}) #,m::mesh,kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    m=p.mesh;
    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    degF=p.degFBoundary;
    cH=projectRecovery(degF[comp[2]],degF[comp[1]],cval,p.massMBoundary[comp[2]],m,kubPoints,kubWeights);
    cHP=projectRecovery(degF[comp[4]],degF[comp[2]],cH,p.massMBoundary[comp[4]],m,kubPoints,kubWeights);
    #=
    n=m.topology.size[3]
    dim=Val{dim}();
    cR=embed(p,comp[2],comp[3],cH,n,dim);
    cEmbed=embed(p,comp[1],comp[3],cval,n,dim);
    cHPEmbed=embed(p,comp[4],comp[3],cHP,n,dim);
    =#
    cR=projectRecovery(degF[comp[3]],degF[comp[2]],cH,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);
    cEmbed=projectRecovery(degF[comp[3]],degF[comp[1]],cval,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);
    cHPEmbed=projectRecovery(degF[comp[3]],degF[comp[4]],cHP,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);

    return cR+(cEmbed-cHPEmbed);
end

function recovery(p::femProblem, dim::Int64, comp::Array{Symbol,1}, cval::Array{Float64,1}, bcomp::Symbol)
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

    #=
    n=m.topology.size[3]
    dim=Val{dim}();
    cR=embed(p,comp[2],comp[3],cH,n,dim);
    cEmbed=embed(p,comp[1],comp[3],cval,n,dim);
    cHPEmbed=embed(p,comp[4],comp[3],cHP,n,dim);
    =#
    cR=projectRecovery(degF[comp[3]],degF[comp[2]],cH,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);
    cEmbed=projectRecovery(degF[comp[3]],degF[comp[1]],cval,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);
    cHPEmbed=projectRecovery(degF[comp[3]],degF[comp[4]],cHP,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);

    return cR+(cEmbed-cHPEmbed);
end
