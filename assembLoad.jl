function assembLoad(degF::degF{1}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=degF.phi;
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB)

    iter=length(phiT);
    for k in 1:m.topology.size[m.topology.dim+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,dJ,m,k,kubPoints,jcoord);
        ft=Array{Float64,2}(undef,sk);
        if sk[1]==1 # <=> dreiecke, muss liste durchlaufen
            for i=1:sk[2]
                xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,i])
                ft[i]=f(xy);
            end
        else # <=> vierecke, muss matrix durchlaufen
            for i=1:sk[1], j=1:sk[2]
                xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
                ft[i,j]=f(xy);
            end
        end


        globalNum=l2g(degF,k);
        for j in 1:iter
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gb[globalNum[j]]+=kubWeights[l,r]*phiT[j][l,r]*ft[l,r]*abs(dJ[l,r]);
                end
            end
        end
    end
    return gb;
end

function assembLoad(degF::degF{2}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=degF.phi;
    sk=size(kubWeights);
    iter=size(phiT,2);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,iter),sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB);

    for k in 1:m.topology.size[m.topology.dim+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,jcoord);
        ft=Array{Array{Float64,2},1}(undef,m.geometry.dim);
        for d in 1:m.geometry.dim
            ft[d]=Array{Float64,2}(undef,sk);
            if sk[1]==1 # <=> dreiecke, muss liste durchlaufen
                for i=1:sk[2]
                    xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,i])
                    ft[d][i]=f[d](xy);
                end
            else # <=> vierecke, muss matrix durchlaufen
                for i=1:sk[1], j=1:sk[2]
                    xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
                    ft[d][i,j]=f[d](xy);
                end
            end
        end
        globalNum=l2g(degF,k);
        for j in 1:iter
            for r in 1:sk[2]
                for l in 1:sk[1]
                    vecdot=0.0
                    for d in 1:m.geometry.dim
                        vecdot+=ft[d][l,r]*jphiT[d,j][l,r]
                    end
                    gb[globalNum[j]]+=kubWeights[l,r]*vecdot*ddJ[l,r]/abs(ddJ[l,r]);
                end
            end
        end
    end
    return gb;
end
