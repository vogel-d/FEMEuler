function projectPressure(degFP::degF{1},massMP::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},
                        degFRT::degF{1},valRT::Array{Float64,1},
                        m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiRT=@views degFRT.phi;
    phiP=@views degFP.phi;
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumRT=Array{Int64,1}(undef,length(phiRT));
    globalNumP=Array{Int64,1}(undef,length(phiP));

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
                    gbh[globalNumP[j]]+=kubWeights[l,r]*cl[l,r]*phiP[j][l,r]*abs(dJ[l,r]);
                end
            end
        end
    end
    return massMP\gbh;
end

function projectPressure(degFP::degF{1},massMP::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},
                        degFRT::degF{1},valRT::Array{Float64,1},
                        m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiRT=@views degFRT.phi;
    phiP=@views degFP.phi;
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumRT=Array{Int64,1}(undef,length(phiRT));
    globalNumP=Array{Int64,1}(undef,length(phiP));

    grav2=0.5*9.81;

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

        @. cl=grav2*cl^2;

        for j in 1:length(phiP)
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumP[j]]+=kubWeights[l,r]*cl[l,r]*phiP[j][l,r]*abs(dJ[l,r]);
                end
            end
        end
    end
    return massMP\gbh;
end
