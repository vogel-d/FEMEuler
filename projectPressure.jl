function projectPressure(degFP::degF{3},massMP::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},
                        degFRT::degF{3},valRT::Array{Float64,1},
                        m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiRT=@views degFRT.phi;
    phiP=@views degFP.phi;
    sk=size(kubWeights);

    #display(phiRT)
    #display(phiP)

    J=initPhi((2,2),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,2,m.meshType);

    globalNumRT=Array{Int64,1}(undef,size(phiRT,3));
    globalNumP=Array{Int64,1}(undef,size(phiP,3));

    cl=zeros(sk);

    gbh=zeros(size(degFP.coordinates,2))
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        l2g!(globalNumRT,degFRT,k);
        l2g!(globalNumP,degFP,k);

        fill!(cl,0.0);
        for i in 1:length(globalNumRT)
            @views @. cl+=valRT[globalNumRT[i]]*phiRT[:,:,i];
        end

        Cpd=1004.0;
        Rd=Cpd-717.0; kappa=Rd/Cpd;
        p0=100000.0;
        Pres1=p0*(Rd/p0)^(1.0/(1.0-kappa));
        Pres2=(1.0/(1.0-kappa));
        cl=Pres1.*cl.^Pres2;

        for j in 1:length(globalNumP)
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumP[j]]+=kubWeights[l,r]*cl[l,r]*phiP[l,r,j]*dJ[l,r];
                end
            end
        end
    end
    return massMP\gbh;
end
