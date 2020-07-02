function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1})
    cS=assembRecovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[1]],cval,p.mesh,p.kubPoints,p.kubWeights)
    cR=recovery(p.degFBoundary[comp[2]],cS,p.stencil,p.recoveryM[(comp[2],comp[3])],p.mesh.topology.size[3])
    return cR
end

function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1},i::Int)
    #cS=assembRecovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[1]],cval,p.mesh,p.kubPoints,p.kubWeights)
    cR=recovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[1]],cval,p.stencil,p.recoveryM[(comp[2],comp[3])],p.mesh,p.kubPoints,p.kubWeights)
    return cR
end

function recovery(degFT::degF{1,:H1}, cval::Array{Float64,1},
                  stencil::Array{Array{Int,1},1}, recoveryM::Array{Any,1}, nf::Int)
    phiT=@views degFT.phi;
    nT=length(phiT)

    globalNumT=Array{Int64,1}(undef,nT);
    lcval=Float64[]
    cR=Float64[];
    for f in 1:nf
        for k in stencil[f]
            k==f ? w=10.0^8 : w=1.0;
            l2g!(globalNumT,degFT,k)
            append!(lcval,w.*cval[globalNumT])
        end
        append!(cR,recoveryM[f]\lcval)
        empty!(lcval);
    end
    return cR;
end

function recovery(degFT::degF{2,:H1div}, degFF::degF{2,:H1div}, cval::Array{Float64,1},
                  stencil::Array{Array{Int,1},1}, recoveryM::Array{Any,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    nT=size(phiT,2)
    nF=size(phiF,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    jphiF=initJacobi((m.geometry.dim,nF),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);
    lcval=Float64[]
    cR=Float64[];
    for f in 1:nf
        z=1;
        for k in stencil[f]
            jacobi!(J,ddJ,jphiT,jphiF,m,k,kubPoints,phiT,phiF,coord);
            k==f ? w=10.0^8 : w=1.0;

            l2g!(globalNumT,degFT,k)
            l2g!(globalNumF,degFF,k);

            for i in 1:nT
                currentval=0.0;
                for j in 1:nF
                    gj=globalNumF[j];
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            vecdot=0.0
                            for d in 1:m.geometry.dim
                                vecdot+=jphiT[d,i][l,r]*jphiF[d,j][l,r]
                            end
                            currentval+=kubWeights[l,r]*ddJ[l,r]^2/abs(ddJ[l,r])*cval[gj]*vecdot;
                        end
                    end
                end
                push!(lcval,w.*currentval)
                z+=1;
            end
        end
        append!(cR,recoveryM[f]\lcval)
        empty!(lcval);
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
