function projectRecovery(degFH::degF{1,:H1},degF::degF{1,:H1},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
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
                    gbh[globalNumH[j]]+=kubWeights[l,k]*cl[l,k]*phiH[j][l,k]*abs(dJ[l,k]);
                end
            end
        end
    end
    return massMH\gbh;
end

function projectRecovery(degFH::degF{2,:H1xH1},degF::degF{2,:H1div},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=size(phiH,2);
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,size(phi,2)),sk);
    jphiH=initJacobi((m.geometry.dim,sph),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    cl=[zeros(sk) for d in 1:m.geometry.dim]

    globalNum=Array{Int64,1}(undef,size(phi,2));
    globalNumH=Array{Int64,1}(undef,size(phiH,2));

    gbh=zeros(degFH.numB)
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphi,m,k,kubPoints,phi,coord);

        l2g!(globalNum,degF,k);
        l2g!(globalNumH,degFH,k);

        for d in 1:m.geometry.dim
            fill!(cl[d], 0.0)
            for i in 1:length(globalNum)
                @. cl[d]+=cval[globalNum[i]]*jphi[d,i];
            end
        end

        for j in 1:sph
            for r in 1:sk[2]
                for l in 1:sk[1]
                    vecdot=0.0
                    for d in 1:m.geometry.dim
                        vecdot+=cl[d][l,r]*phiH[d,j][l,r];
                    end
                    gbh[globalNumH[j]]+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*vecdot;
                    #piola: abs(dJ)*ddJ=ddJ/abs(ddJ)
                end
            end
        end
    end
    return massMH\gbh;
end

function projectRecovery(degFH::degF{2,:H1div},degF::degF{2,:H1xH1},cval::Array{Float64,1},massMH::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64},m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiH=@views degFH.phi;
    sph=size(phiH,2);
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,size(phi,2)),sk);
    jphiH=initJacobi((m.geometry.dim,sph),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    cl=[zeros(sk) for d in 1:m.geometry.dim]

    globalNum=Array{Int64,1}(undef,size(phi,2));
    globalNumH=Array{Int64,1}(undef,size(phiH,2));

    gbh=zeros(degFH.numB)
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphiH,m,k,kubPoints,phiH,coord);

        l2g!(globalNum,degF,k);
        l2g!(globalNumH,degFH,k);

        for d in 1:m.geometry.dim
            fill!(cl[d], 0.0)
            for i in 1:length(globalNum)
                @. cl[d]+=cval[globalNum[i]]*phi[d,i];
            end
        end

        for j in 1:sph
            for r in 1:sk[2]
                for l in 1:sk[1]
                    vecdot=0.0
                    for d in 1:m.geometry.dim
                        vecdot+=cl[d][l,r]*jphiH[d,j][l,r];
                    end
                    gbh[globalNumH[j]]+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*vecdot;
                    #piola: abs(dJ)*ddJ=ddJ/abs(ddJ)
                end
            end
        end
    end
    return massMH\gbh;
end
