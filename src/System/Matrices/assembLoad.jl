function assembLoad(degF::degF{1,:H1}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=degF.phi;
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB)

    ft=zeros(sk);

    iter=length(phiT);
    for k in 1:m.topology.size[m.topology.dim+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,dJ,m,k,kubPoints,jcoord);
        fill!(ft,0.0);
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

function assembLoad(degF::degF{2,:H1div}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=degF.phi;
    sk=size(kubWeights);
    iter=size(phiT,2);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,iter),sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB);

    ft=[zeros(sk) for d in 1:m.geometry.dim]

    for k in 1:m.topology.size[m.topology.dim+1]
        coord=@views m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,jcoord);
        for d in 1:m.geometry.dim
            fill!(ft[d],0.0);
            if sk[1]==1 # <=> dreiecke, muss liste durchlaufen
                for i=1:sk[2]
                    xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,i])
                    ft[d][i]=f(xy)[d];
                end
            else # <=> vierecke, muss matrix durchlaufen
                for i=1:sk[1], j=1:sk[2]
                    xy=transformation(m,coord,kubPoints[1,i],kubPoints[2,j])
                    ft[d][i,j]=f(xy)[d];
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


function assembLoad!(F::Array{Float64,1},
                     degFRT::degF{2,:H1div}, f::Array{Float64,1},
                     degFVecDG::degF{2,:H1xH1},
                     m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})

    phi=degFRT.phi;
    psi=degFVecDG.phi;
    nphi=size(phi,2);
    npsi=size(psi,2);

    gverticesRT=Array{Int64,1}(undef,nphi)
    gverticesVecDG=Array{Int64,1}(undef,npsi)

    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,npsi),sk);

    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vecdot=0.0;
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphi,m,k,kubPoints,phi,coord);

        l2g!(gverticesRT,degFRT,k);
        l2g!(gverticesVecDG,degFVecDG,k);

        for i in 1:npsi
            for j in 1:nphi
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        vecdot=0.0;
                        for d in 1:m.geometry.dim
                            vecdot+=psi[d,i][l,r]*jphi[d,j][l,r]
                        end
                        F[gverticesVecDG[i]]+=kubWeights[l,r]*f[gverticesRT[j]]*(ddJ[l,r]/abs(ddJ[l,r]))*vecdot;
                    end
                end
            end
        end
    end
end
