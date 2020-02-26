function projectRecovery(degFH::degF{1},degF::degF{1},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=length(phiH)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    cl=zeros(sk);

    globalNum=Array{Int64,1}(undef,length(phi));
    globalNumH=Array{Int64,1}(undef,length(phiH));

    gbh=zeros(degFH.numB)
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        l2g!(globalNum,degF,k);
        l2g!(globalNumH,degFH,k);

        fill!(cl,0.0);
        for i in 1:length(globalNum)
            @. cl+=cval[globalNum[i]]*phi[i];
        end

        for j in 1:sph
            for k in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumH[j]]+=kubWeights[l,k]*cl[l,k]*phiH[j][l,k]*dJ[l,k];
                end
            end
        end
    end
    return massMH\gbh;
end

function projectRecovery(degFH::degF{2},degF::degF{2},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=size(phiH,2);
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,size(phi,2)),sk);
    jphiH=initJacobi((m.geometry.dim,sph),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    cl=Array{Array{Float64,2},1}(undef,m.geometry.dim);
    fill!(cl, zeros(sk))

    globalNum=Array{Int64,1}(undef,size(phi,2));
    globalNumH=Array{Int64,1}(undef,size(phiH,2));

    gbh=zeros(degFH.numB)
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphi,jphiH,m,k,kubPoints, phi, phiH,coord);

        l2g!(globalNum,degF,k);
        l2g!(globalNumH,degFH,k);

        fill!(cl, zeros(sk))
        for i in 1:length(globalNum)
            for d in 1:m.geometry.dim
                @. cl[d]+=cval[globalNum[i]]*jphi[d,i];
            end
        end

        for j in 1:sph
            for r in 1:sk[2]
                for l in 1:sk[1]
                    vecdot=0.0
                    for d in 1:m.geometry.dim
                        vecdot+=cl[d][l,r]*jphiH[d,j][l,r];
                    end
                    gbh[globalNumH[j]]+=kubWeights[l,r]*ddJ[l,r]*vecdot;
                end
            end
        end
    end
    return massMH\gbh;
end
