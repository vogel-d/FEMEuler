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
                ft[i]=f(xy[1],xy[2]);
            end
        else # <=> vierecke, muss matrix durchlaufen
            for i=1:sk[1], j=1:sk[2]
                xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
                ft[i,j]=f(xy[1],xy[2]);
            end
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
    iter=size(phiT,2);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,iter),sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB);

    for k in 1:m.topology.size[m.topology.dim+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,jcoord);
<<<<<<< HEAD
        ft=Array{Array{Float64,2},1}(undef,m.geometry.dim);
        for d in 1:m.geometry.dim
            ft[d]=Array{Float64,2}(undef,sk);
            for i=1:sk[1], j=1:sk[2]
                xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
                ft[d][i,j]=f[1](xy[1],xy[2]);
=======
        ft1=Array{Float64,2}(undef,sk);
        ft2=Array{Float64,2}(undef,sk);
        if sk[1]==1 # <=> dreiecke, muss liste durchlaufen
            for i=1:sk[2]
                xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,i])
                ft1[i]=f[1](xy[1],xy[2]);
                ft2[i]=f[2](xy[1],xy[2]);
            end
        else # <=> vierecke, muss matrix durchlaufen
            for i=1:sk[1], j=1:sk[2]
                xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
                ft1[i,j]=f[1](xy[1],xy[2]);
                ft2[i,j]=f[2](xy[1],xy[2]);
>>>>>>> 9769bfbcb07b5b2b7e36bbd8673a21c1644bd528
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
                    gb[globalNum[j]]+=kubWeights[l,r]*vecdot;
                end
            end
        end
    end
    return gb;
end
