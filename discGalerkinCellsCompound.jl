function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, divphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=size(kubWeights);

    nSubCells=compoundData.nSubCells;
    w=zeros(sk);
    gradw1=zeros(sk);
    gradw2=zeros(sk);
    mt=m.meshType;

    nSubPhiT=length(phiT);
    nSubPhiF=size(phiF,2);
    nSubPhiW=length(phiW);
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    #edge integrating for compound RT0
    nquadPhiF=compoundData.nquadPhi[degFF.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiF_edge=initJacobi((m.geometry.dim,size(phiF,2)),sq);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhiT, subcoord, degFT, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF, subcoord, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW, subcoord, degFW, m, J, dJ, phiW, kubPoints, kubWeights, compoundData);

        for subCell in 1:nSubCells
            fill!(w,0.0);
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[i];
                end
            end

            fill!(gradw1,0.0);
            fill!(gradw2,0.0);
            for j in 1:nCompoundPhiW
                for subj in 1:nSubPhiW
                    @. gradw1+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,j]
                    @. gradw2+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,j];
                end
            end

            jacobi!(J,dJ,kubPoints,subcoord[subCell],mt);

            l2g!(globalNumF,degFF,k);
            l2g!(globalNumT,degFT,k);

            for i in 1:nCompoundPhiT
                gi=globalNumT[i];
                z=0.0;
                for j in 1:nCompoundPhiF
                    gj=globalNumF[j];
                    for subi in 1:nSubPhiT
                        for subj in 1:nSubPhiF
                            for r in 1:size(kubWeights,2)
                                for l in 1:size(kubWeights,1)
                                    z+= assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                        fval[gj]*kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*phiT[subi][l,r]*(w[l,r]*divphiF[subj][l,r]+(gradw1[l,r]*phiF[1,subj][l,r]+gradw2[l,r]*phiF[2,subj][l,r]));
                                end
                            end
                        end
                    end
                end
                M[gi]-=z;
            end
        end
    end

    return nothing;
end

function discGalerkinCells!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, divphiF::Array{Array{Float64,2},1}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=size(kubWeights);

    w=zeros(sk);
    gradw1=zeros(sk);
    gradw2=zeros(sk);
    mt=m.meshType;

    nSubCells=compoundData.nSubCells;
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    nSubPhiT=length(phiT);
    nSubPhiF=size(phiF,2);
    nSubPhiW=length(phiW);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    #edge integrating for compound RT0
    nquadPhiF=compoundData.nquadPhi[degFF.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiF_edge=initJacobi((m.geometry.dim,size(phiF,2)),sq);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    lM=zeros(nCompoundPhiT,nCompoundPhiF);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhiT, subcoord, degFT, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF, subcoord, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW, subcoord, degFW, m, J, dJ, phiW, kubPoints, kubWeights, compoundData);

        for subCell in 1:nSubCells

            fill!(w,0.0);
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[i];
                end
            end

            fill!(gradw1,0.0);
            fill!(gradw2,0.0);
            for j in 1:nCompoundPhiW
                for subj in 1:nSubPhiW
                    @. gradw1+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,j]
                    @. gradw2+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,j];
                end
            end

            jacobi!(J,dJ,kubPoints,subcoord[subCell],mt);

            fill!(lM,0.0);
            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    for subj in 1:nSubPhiF
                        for subi in 1:nSubPhiT
                            for r in 1:sk[2]
                                for l in 1:sk[1]
                                    lM[i,j]+=assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                             kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*phiT[subi][l,r]*(w[l,r]*divphiF[subj][l,r]+(gradw1[l,r]*phiF[1,subj][l,r]+gradw2[l,r]*phiF[2,subj][l,r]));
                                end
                            end
                        end
                    end
                end
            end
            l2g!(globalNumF,degFF,k);
            l2g!(globalNumT,degFT,k);

            for i in 1:length(globalNumT)
                gi=globalNumT[i];
                for j in 1:length(globalNumF)
                    push!(rows,gi);
                    push!(cols,globalNumF[j]);
                    push!(vals,-lM[i,j]);
                end
            end
        end
    end

    return nothing;
end


function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, divphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            #degFW::degF{2,S} where S,phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=size(kubWeights);

    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi(size(phiT),sk)

    nSubCells=compoundData.nSubCells;
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    w1=zeros(sk);
    w2=zeros(sk);
    gradw11=zeros(sk);
    gradw12=zeros(sk);
    gradw21=zeros(sk);
    gradw22=zeros(sk);
    mt=m.meshType;

    nSubPhiT=size(phiT,2);
    nSubPhiF=size(phiF,2);
    nSubPhiW=size(phiW,2);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    #edge integrating for compound RT0
    divphiT=degFT.divphi;
    divphiW=degFW.divphi;
    nquadPhiT=compoundData.nquadPhi[degFT.femType];
    nquadPhiF=compoundData.nquadPhi[degFF.femType];
    nquadPhiW=compoundData.nquadPhi[degFW.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiT_edge=initJacobi((m.geometry.dim,size(phiT,2)),sq);
    jphiF_edge=initJacobi((m.geometry.dim,size(phiF,2)),sq);
    jphiW_edge=initJacobi((m.geometry.dim,size(phiW,2)),sq);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhiT, subcoord, m, divphiT, J_edge, ddJ_edge, jphiT_edge, nquadPhiT, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiF, subcoord, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW, subcoord, m, divphiW, J_edge, ddJ_edge, jphiW_edge, nquadPhiW, nquadPoints, quadWeights, compoundData);

        l2g!(globalNumW,degFW,k);
        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for subCell in 1:nSubCells
            fill!(w1,0.0);
            fill!(w2,0.0);
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1+=assembledPhiW[i][subi,subCell]*
                           wval[globalNumW[i]]*phiW[1,subi];
                    @. w2+=assembledPhiW[i][subi,subCell]*
                           wval[globalNumW[i]]*phiW[2,subi];
                end
            end

            fill!(gradw11,0.0);
            fill!(gradw12,0.0);
            fill!(gradw21,0.0);
            fill!(gradw22,0.0);
            for j in 1:nCompoundPhiW
                for subj in 1:nSubPhiW
                    @. gradw11+=assembledPhiW[j][subj,subCell]*
                                wval[globalNumW[j]]*gradphiW[1,1+2*(subj-1)]
                    @. gradw12+=assembledPhiW[j][subj,subCell]*
                                wval[globalNumW[j]]*gradphiW[1,2+2*(subj-1)];
                    @. gradw21+=assembledPhiW[j][subj,subCell]*
                                wval[globalNumW[j]]*gradphiW[2,1+2*(subj-1)];
                    @. gradw22+=assembledPhiW[j][subj,subCell]*
                                wval[globalNumW[j]]*gradphiW[2,2+2*(subj-1)];

                end
            end

            jacobi!(ddJ,jphiT,kubPoints,phiT,subcoord[subCell],mt);

            for i in 1:nCompoundPhiT
                gi=globalNumT[i];
                z=0.0;
                for j in 1:nCompoundPhiF
                    gj=globalNumF[j];
                    for subi in 1:nSubPhiT
                        for subj in 1:nSubPhiF
                            for r in 1:size(kubWeights,2)
                                for l in 1:size(kubWeights,1)
                                    z+=assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                       fval[gj]*kubWeights[l,r]*(ddJ[l,r]^3/abs(ddJ[l,r]))*(jphiT[1,subi][l,r]*(divphiF[subj][l,r]*w1[l,r]+gradw11[l,r]*phiF[1,subj][l,r]+gradw12[l,r]*phiF[2,subj][l,r])+jphiT[2,subi][l,r]*(divphiF[subj][l,r]*w2[l,r]+gradw21[l,r]*phiF[1,subj][l,r]+gradw22[l,r]*phiF[2,subj][l,r]));
                                end
                            end
                        end
                    end
                end
                M[gi]-=z;
            end
        end
    end
    return nothing;
end


function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, divphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1xH1},phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)
    @error("VecDG1 not assembled for compound elements")
    sk=size(kubWeights);
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi(size(phiT),sk);
    jphiF=initJacobi(size(phiF),sk);

    nSubCells=compoundData.nSubCells;
    w_clean=zeros(sk);
    w1=Array{Array{Float64,2},1}(undef,nSubCells);
    w2=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw11=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw12=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw21=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw22=Array{Array{Float64,2},1}(undef,nSubCells);
    mt=m.meshType;

    nSubPhiW=length(phiW);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    for k in 1:m.topology.size[3]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhi, subcoord, degF, m, J, dJ, phi, kubPoints, kubWeights, compoundData);

        fill!(w1,w_clean);
        fill!(w2,w_clean);
        for subCell in 1:nSubCells
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[1,i];
                    @. w2+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[2,i];
                end
            end
        end

        fill!(gradw11,w_clean);
        fill!(gradw12,w_clean);
        fill!(gradw21,w_clean);
        fill!(gradw22,w_clean);
        for subCell in 1:nSubCells
            for j in 1:nCompoundPhiW
                for subj in 1:nSubPhiW
                    @. gradw11+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,1+2*(j-1)]
                    @. gradw12+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,2+2*(j-1)];
                    @. gradw21+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,1+2*(j-1)];
                    @. gradw22+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,2+2*(j-1)];

                end
            end
        end

        l2g!(globalNumW,degFW,k);
        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for subCell in 1:nSubCells
            jacobi!(J,ddJ,jphiT,kubPoints,phiT,subcoord[subCell],mt);
            for i in 1:nCompoundPhiT
                gi=globalNumT[i];
                z=0.0;
                for j in 1:nCompoundPhiF
                    gj=globalNumF[j];
                    for subi in nSubPhiT
                        for subj in 1:nSubPhiF
                            for r in 1:size(kubWeights,2)
                                for l in 1:size(kubWeights,1)
                                    z+=assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                       fval[gj]*kubWeights[l,r]*abs(ddJ[l,r])*(jphiT[1,i][l,r]*(divphiF[subj][l,r]*w1[l,r]+gradw11[l,r]*phiF[1,subj][l,r]+gradw12[l,r]*phiF[2,subj][l,r])+jphiT[2,subi][l,r]*(divphiF[subj][l,r]*w1[l,r]+gradw21[l,r]*phiF[1,subj][l,r]+gradw22[l,r]*phiF[2,subj][l,r]));
                                end
                            end
                        end
                    end
                end
                M[gi]-=z;
            end
        end
    end
    return nothing;
end
