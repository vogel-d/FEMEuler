function projectRecovery(degFH::degF{1},degF::degF{1},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=length(phiH)
    sk=size(kubWeights);

    J=Array{Array{Float64,2},2}(undef,2,2);
    dJ=Array{Float64,2}(undef,sk);

    gbh=zeros(size(degFH.coordinates,2))
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints);

        globalNum=@views l2g(degF,k);
        globalNumH=@views l2g(degFH,k);

        cl=zeros(sk);
        for i in 1:length(globalNum)
            cl+=cval[globalNum[i]]*phi[i];
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

    J=Array{Array{Float64,2},2}(undef,2,2);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initPhi(size(phi),sk);
    jphiH=initPhi(size(phiH),sk);

    gbh=zeros(size(degFH.coordinates,2))
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphi,jphiH,m,k,kubPoints, phi, phiH);

        globalNum=@views l2g(degF,k);
        globalNumH=@views l2g(degFH,k);

        cl1=zeros(sk);
        cl2=zeros(sk);
        for i in 1:length(globalNum)
            cl1+=cval[globalNum[i]]*jphi[1,i];
            cl2+=cval[globalNum[i]]*jphi[2,i];
        end

        for j in 1:sph
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumH[j]]+=kubWeights[l,r]*ddJ[l,r]*(cl1[l,r]*jphiH[1,j][l,r]+cl2[l,r]*jphiH[2,j][l,r]);
                end
            end
        end
    end
    return massMH\gbh;
end
