function assembRecovery(degFT::degF{1,:H1},degFF::degF{1,:H1},valF::Array{Float64,1},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    nT=length(phiT)
    nF=length(phiF)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);

    cL=zeros(degFT.numB)
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        l2g!(globalNumT,degFT,k);
        l2g!(globalNumF,degFF,k);

        for i in 1:nT
            gi=globalNumT[i];
            z=0.0;
            for j in 1:nF
                gj=globalNumF[j];
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        z+=phiT[i][l,r]*kubWeights[l,r]*abs(dJ[l,r])*valF[gj]*phiF[j][l,r];
                    end
                end
            end
            cL[gi] += z;
        end
    end
    return cL;
end
#=
function assembRecovery(degFT::degF{2,:H1div},degFF::degF{2,:H1div},valF::Array{Float64,1},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
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

    cL=zeros(degFT.numB)
    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiT,jphiF,m,k,kubPoints,phiT,phiF,coord);

        l2g!(globalNumT,degFT,k);
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
                        z+=kubWeights[l,r]*ddJ[l,r]^2/abs(ddJ[l,r])*valF[gj]*vecdot;
                    end
                end
            end
            cL[gi] += z;
        end
    end
    return cL;
end

function assembRecovery(degFT::degF{2,:H1xH1},degFF::degF{2,:H1xH1},valF::Array{Float64,1},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    nT=size(phiT,2)
    nF=size(phiF,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);

    cL=zeros(degFT.numB)
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        l2g!(globalNumT,degFT,k);
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
                            vecdot+=phiT[d,i][l,r]*phiF[d,j][l,r]
                        end
                        z+=kubWeights[l,r]*abs(dJ[l,r])*valF[gj]*vecdot;
                    end
                end
            end
            cL[gi] += z;
        end
    end
    return cL;
end


function assembRecovery(degFT::degF{2,:H1div},degFF::degF{2,:H1xH1},valF::Array{Float64,1},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    nT=size(phiT,2)
    nF=size(phiF,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);

    cL=zeros(degFT.numB)
    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,coord);

        l2g!(globalNumT,degFT,k);
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
                            vecdot+=jphiT[d,i][l,r]*phiF[d,j][l,r]
                        end
                        z+=kubWeights[l,r]*ddJ[l,r]/abs(ddJ[l,r])*valF[gj]*vecdot;
                    end
                end
            end
            cL[gi] += z;
        end
    end
    return cL;
end

function assembRecovery(degFT::degF{2,:H1xH1},degFF::degF{2,:H1div},valF::Array{Float64,1},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    nT=size(phiT,2)
    nF=size(phiF,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiF=initJacobi((m.geometry.dim,nF),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nT);
    globalNumF=Array{Int64,1}(undef,nF);

    cL=zeros(degFT.numB)
    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiF,m,k,kubPoints,phiF,coord);

        l2g!(globalNumT,degFT,k);
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
                            vecdot+=phiT[d,i][l,r]*jphiF[d,j][l,r]
                        end
                        z+=kubWeights[l,r]*ddJ[l,r]/abs(ddJ[l,r])*valF[gj]*vecdot;
                    end
                end
            end
            cL[gi] += z;
        end
    end
    return cL;
end
=#
