function projectPressure(degFP::degF{1},massMP::SuiteSparse.UMFPACK.UmfpackLU{AbstractFloat,Int},
                        degFRT::degF{1},valRT::Array{AbstractFloat,1},
                        m::mesh, kubPoints::Array{AbstractFloat,2}, kubWeights::Array{AbstractFloat,2})
    phiRT=@views degFRT.phi;
    phiP=@views degFP.phi;
    sk=size(kubWeights);

    J=initJacobi((2,2),sk);
    dJ=Array{AbstractFloat,2}(undef,sk);
    coord=Array{AbstractFloat,2}(undef,2,m.meshType);

    globalNumRT=Array{Int,1}(undef,length(phiRT));
    globalNumP=Array{Int,1}(undef,length(phiP));

    Cpd=1004.0;
    Rd=Cpd-717.0; kappa=Rd/Cpd;
    p0=100000.0;
    Pres1=p0*(Rd/p0)^(1.0/(1.0-kappa));
    Pres2=(1.0/(1.0-kappa));

    cl=zeros(sk);

    gbh=zeros(degFP.numB)
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        l2g!(globalNumRT,degFRT,k);
        l2g!(globalNumP,degFP,k);

        fill!(cl,0.0);
        for i in 1:length(globalNumRT)
            @. cl+=valRT[globalNumRT[i]]*phiRT[i];
        end

        @. cl=Pres1*cl^Pres2;

        for j in 1:length(phiP)
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumP[j]]+=kubWeights[l,r]*cl[l,r]*phiP[j][l,r]*dJ[l,r];
                end
            end
        end
    end
    return massMP\gbh;
end
