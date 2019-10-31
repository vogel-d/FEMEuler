function assembLoad(degF::degF{1}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=degF.phi;
    sk=size(kubWeights);

    J=initPhi((2,2),sk);
    dJ=Array{Float64,2}(undef,sk);
    jcoord=Array{Float64,2}(undef,2,m.meshType);
    gb=zeros(size(degF.coordinates,2))

    iter=length(phiT);
    for k in 1:m.topology.size[m.topology.D+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,dJ,m,k,kubPoints,jcoord);
        ft=Array{Float64,2}(undef,sk);
        for i=1:sk[1], j=1:sk[2]
            xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
            ft[i,j]=f(xy[1],xy[2]);
        end

        globalNum=l2g(degF,k);
        for j in 1:iter
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gb[globalNum[j]]+=kubWeights[l,r]*phiT[j][l,r]*ft[l,r]*dJ[l,r];
                end
            end
        end
    end
    return gb;
end

function assembLoad(degF::degF{2}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=degF.phi;
    sk=size(kubWeights);

    J=initPhi((2,2),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initPhi(size(phiT),sk);
    jcoord=Array{Float64,2}(undef,2,m.meshType);
    gb=zeros(size(degF.coordinates,2));

    iter=size(phiT,2);
    for k in 1:m.topology.size[m.topology.D+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,jcoord);
        ft1=Array{Float64,2}(undef,sk);
        ft2=Array{Float64,2}(undef,sk);
        for i=1:sk[1], j=1:sk[2]
            xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
            ft1[i,j]=f[1](xy[1],xy[2]);
            ft2[i,j]=f[2](xy[1],xy[2]);
        end

        globalNum=l2g(degF,k);
        for j in 1:iter
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gb[globalNum[j]]+=kubWeights[l,r]*(ft1[l,r]*jphiT[1,j][l,r]+ft2[l,r]*jphiT[2,j][l,r]);
                end
            end
        end
    end
    return gb;
end
