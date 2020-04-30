function assembLoadCompound(degF::degF{1,:H1}, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, compoundData::compoundData)
    phi=degF.phi;
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB)

    nSubCells=compoundData.nSubCells;
    nCompoundPhi=compoundData.nCompoundPhi[degF.femType];
    assembledPhi=compoundData.assembledPhi[degF.femType];
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);

    ft=zeros(sk);

    nPhiSubElement=length(phi);
    for k in 1:m.topology.size[m.topology.dim+1]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        getSubCells!(subcoord, coord, compoundData);
        globalNum=l2g(degF,k);
        assemblePhi!(assembledPhi, subcoord, degF, m, J, dJ, phi, kubPoints, kubWeights, compoundData);
        for subCell in 1:nSubCells
            jacobi!(J,dJ,kubPoints,subcoord[subCell],m.meshType);
            for j in 1:nCompoundPhi
                for d in 1:m.geometry.dim
                    fill!(ft,0.0);
                    if sk[1]==1 # <=> dreiecke, muss liste durchlaufen
                        for i=1:sk[2]
                            xy=transformation(m,subcoord[subCell],kubPoints[1,i],kubPoints[2,i])
                            ft[i]=f(xy);
                        end
                    else # <=> vierecke, muss matrix durchlaufen
                        for i=1:sk[1], j=1:sk[2]
                            xy=transformation(m,subcoord[subCell],kubPoints[1,i],kubPoints[2,j])
                            ft[i,j]=f(xy);
                        end
                    end
                end

                for subj in 1:nPhiSubElement
                    if assembledPhi[j][subj,subCell]!=0
                        for r in 1:sk[2]
                            for l in 1:sk[1]
                                gb[globalNum[j]]+=assembledPhi[j][subj,subCell]*
                                                  kubWeights[l,r]*phi[subj][l,r]*ft[l,r]*abs(dJ[l,r]);
                            end
                        end
                    end
                end
            end
        end
    end
    return gb;
end

function assembLoadCompound(degF::degF{2,S} where S, f, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, compoundData::compoundData)
    phi=degF.phi;
    sk=size(kubWeights);
    nPhiSubElement=size(phi,2);

    coord=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,iter),sk);
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    gb=zeros(degF.numB);
    sq=length(quadWeights)
    J1=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ1=Array{Float64,1}(undef,sq);
    jphi1=initJacobi((m.geometry.dim,size(phi,2)),sq);

    nSubCells=compoundData.nSubCells;
    assembledPhi=compoundData.assembledPhi[degF.femType];
    nCompoundPhi=compoundData.nCompoundPhi[degF.femType];
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);

    nquadPhi=compoundData.nquadPhi[degF.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;

    ft=[zeros(sk) for d in 1:m.geometry.dim]

    for k in 1:m.topology.size[m.topology.dim+1]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        getSubCells!(subcoord, coord, compoundData);
        globalNum=l2g(degF,k);
        assemblePhi!(assembledPhi, subcoord, m, divphi, J1, ddJ1, jphi1, nquadPhi, nquadPoints, quadWeights, compoundData);
        for subCell in 1:nSubCells
            jacobi!(J,ddJ,jphi,kubPoints,phi,subcoord[subCell],m.meshType);
            for j in 1:nCompoundPhi
                for d in 1:m.geometry.dim
                    fill!(ft[d],0.0);
                    if sk[1]==1 # <=> dreiecke, muss liste durchlaufen
                        for i=1:sk[2]
                            xy=transformation(m,subcoord[subCell],kubPoints[1,i],kubPoints[2,i])
                            ft[d][i]=f[d](xy);
                        end
                    else # <=> vierecke, muss matrix durchlaufen
                        for i=1:sk[1], j=1:sk[2]
                            xy=transformation(m,subcoord[subCell],kubPoints[1,i],kubPoints[2,j])
                            ft[d][i,j]=f[d](xy);
                        end
                    end
                end

                for subj in 1:nPhiSubElement
                    if assembledPhi[j][subj,subCell]!=0
                        for r in 1:sk[2]
                            for l in 1:sk[1]
                                vecdot=0.0
                                for d in 1:m.geometry.dim
                                    vecdot+=ft[d][l,r]*jphi[d,subj][l,r]
                                end
                                gb[globalNum[j]]+=assembledPhi[j][subj,subCell]*
                                                  kubWeights[l,r]*vecdot*ddJ[l,r]/abs(ddJ[l,r]);
                            end
                        end
                    end
                end
            end
        end
    end
    return gb;
end


function assembLoadCompound!(F::Array{Float64,1},
                     degFRT::degF{2,:H1div}, f::Array{Float64,1},
                     degFVecDG::degF{2,:H1xH1},
                     m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, compoundData::compoundData)

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
