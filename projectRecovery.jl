function projectRecovery(degFH::degF{3},degF::degF{3},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=size(phiH,3)
    sk=size(kubWeights);

    J=initPhi((2,2),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,2,m.meshType);

    cl=zeros(sk);

    globalNum=Array{Int64,1}(undef,size(phi,3));
    globalNumH=Array{Int64,1}(undef,size(phiH,3));

    gbh=zeros(size(degFH.coordinates,2))
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        l2g!(globalNum,degF,k);
        l2g!(globalNumH,degFH,k);

        fill!(cl,0.0);
        for i in 1:length(globalNum)
            #@views @. cl+=cval[globalNum[i]]*phi[:,:,i];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    cl[l,r]+=cval[globalNum[i]]*phi[l,r,i];
                end
            end
        end

        for j in 1:sph
            for k in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumH[j]]+=kubWeights[l,k]*cl[l,k]*phiH[l,k,j]*dJ[l,k];
                end
            end
        end
    end
    return massMH\gbh;
end

function projectRecovery(degFH::degF{4},degF::degF{4},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=size(phiH,4);
    sk=size(kubWeights);

    J=initPhi((2,2),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=Array{Float64,4}(undef,size(phi,1),sk[1],sk[2],size(phi,4));
    jphiH=Array{Float64,4}(undef,size(phiH,1),sk[1],sk[2],size(phiH,4));
    coord=Array{Float64,2}(undef,2,m.meshType);

    cl1=zeros(sk);
    cl2=zeros(sk);

    globalNum=Array{Int64,1}(undef,size(phi,4));
    globalNumH=Array{Int64,1}(undef,size(phiH,4));

    gbh=zeros(size(degFH.coordinates,2))
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphi,jphiH,m,k,kubPoints, phi, phiH,coord);

        l2g!(globalNum,degF,k);
        l2g!(globalNumH,degFH,k);

        fill!(cl1,0.0);
        fill!(cl2,0.0);
        for i in 1:length(globalNum)
            #@views @. cl1+=cval[globalNum[i]]*jphi[1,:,:,i];
            #@views @. cl2+=cval[globalNum[i]]*jphi[2,:,:,i];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    cl1[l,r]+=cval[globalNum[i]]*jphi[1,l,r,i];
                    cl2[l,r]+=cval[globalNum[i]]*jphi[2,l,r,i];
                end
            end
        end

        for j in 1:sph
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumH[j]]+=kubWeights[l,r]*ddJ[l,r]*(cl1[l,r]*jphiH[1,l,r,j]+cl2[l,r]*jphiH[2,l,r,j]);
                end
            end
        end
    end
    return massMH\gbh;
end
