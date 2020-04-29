function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=size(kubWeights);

    nSubCells=compoundData.nSubCells;
    w_clean=zeros(sk);
    w=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw1=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw2=Array{Array{Float64,2},1}(undef,nSubCells);
    mt=m.meshType;

    nSubPhiW=length(phiW);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhi, subcoord, degF, m, J, dJ, phi, kubPoints, kubWeights, compoundData);

        fill!(w,w_clean);
        for subCell in 1:nSubCells
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w[subCell]+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[i];
                end
            end
        end

        fill!(gradw1,w_clean);
        fill!(gradw2,w_clean);
        for subCell in 1:nSubCells
            for j in 1:nCompoundPhiW
                for subj in 1:nSubPhiW
                    @. gradw1[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,j]
                    @. gradw2[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,j];
                end
            end
        end

        for subCell in 1:nSubCells
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
                                        fval[gj]*kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*phiT[subi][l,r]*(w[subCell][l,r]*dphiF[subj][l,r]+(gradw1[subCell][l,r]*phiF[1,subj][l,r]+gradw2[subCell][l,r]*phiF[2,subj][l,r]));
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
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=size(kubWeights);

    w_clean=zeros(sk);
    w=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw1=Array{Array{Float64,2},1}(undef,nSubCells);
    gradw2=Array{Array{Float64,2},1}(undef,nSubCells);
    mt=m.meshType;

    nSubCells=compoundData.nSubCells;
    nSubPhiW=length(phiW);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    lM=zeros(nCompoundPhiT,nCompoundPhiF);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhi, subcoord, degF, m, J, dJ, phi, kubPoints, kubWeights, compoundData);

        fill!(w,w_clean);
        for subCell in 1:nSubCells
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w[subCell]+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[i];
                end
            end
        end

        fill!(gradw1,w_clean);
        fill!(gradw2,w_clean);
        for subCell in 1:nSubCells
            for j in 1:nCompoundPhiW
                for subj in 1:nSubPhiW
                    @. gradw1[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,j]
                    @. gradw2[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,j];
                end
            end
        end

        for subCell in 1:nSubCells
            jacobi!(J,dJ,kubPoints,subcoord[subCell],mt);

            fill!(lM,0.0);
            for j in 1:nCompoundPhiF
                for i in 1:nComopundPhiT
                    for subj in 1:nSubPhiF
                        for subi in 1:nSubPhiT
                            for r in 1:sk[2]
                                for l in 1:sk[1]
                                    lM[i,j]+=assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                             kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*phiT[subi][l,r]*(w[subCell][l,r]*dphiF[subj][l,r]+(gradw1[subCell][l,r]*phiF[1,subj][l,r]+gradw2[subCell][l,r]*phiF[2,subj][l,r]));
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
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            #degFW::degF{2,S} where S,phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=size(kubWeights);

    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi(size(phiT),sk)

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
        l2g!(globalNumW,degFW,k);
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, compoundData);
        assemblePhi!(assembledPhi, subcoord, degF, m, J, dJ, phi, kubPoints, kubWeights, compoundData);

        fill!(w1,w_clean);
        fill!(w2,w_clean);
        for subCell in 1:nSubCells
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1[subCell]+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[1,i];
                    @. w2[subCell]+=assembledPhiW[i][subi,subCell]*
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
                    @. gradw11[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,1+2*(j-1)]
                    @. gradw12[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,2+2*(j-1)];
                    @. gradw21[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,1+2*(j-1)];
                    @. gradw22[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,2+2*(j-1)];

                end
            end
        end

        l2g!(globalNumW,degFW,k);
        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for subCell in 1:nSubCells
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
                                       fval[gj]*kubWeights[l,r]*(ddJ[l,r]^3/abs(ddJ[l,r]))*(jphiT[1,subi][l,r]*(dphiF[subj][l,r]*w1[subCell][l,r]+gradw11[subCell][l,r]*phiF[1,subj][l,r]+gradw12[subCell][l,r]*phiF[2,subj][l,r])+jphiT[2,subi][l,r]*(dphiF[subj][l,r]*w2[subCell][l,r]+gradw21[subCell][l,r]*phiF[1,subj][l,r]+gradw22[subCell][l,r]*phiF[2,subj][l,r]));
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
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1xH1},phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2}
                            compoundData::compoundData)

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
                    @. w1[subCell]+=assembledPhiW[i][subi,subCell]*
                                   wval[globalNumW[i]]*phiW[1,i];
                    @. w2[subCell]+=assembledPhiW[i][subi,subCell]*
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
                    @. gradw11[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,1+2*(j-1)]
                    @. gradw12[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[1,2+2*(j-1)];
                    @. gradw21[subCell]+=assembledPhiW[j][subj,subCell]*
                                        wval[globalNumW[j]]*gradphiW[2,1+2*(j-1)];
                    @. gradw22[subCell]+=assembledPhiW[j][subj,subCell]*
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
                                       fval[gj]*kubWeights[l,r]*abs(ddJ[l,r])*(jphiT[1,i][l,r]*(dphiF[subj][l,r]*w1[subCell][l,r]+gradw11[subCell][l,r]*phiF[1,subj][l,r]+gradw12[subCell][l,r]*phiF[2,subj][l,r])+jphiT[2,subi][l,r]*(dphiF[subj][l,r]*w2[subCell][l,r]+gradw21[subCell][l,r]*phiF[1,subj][l,r]+gradw22[subCell][l,r]*phiF[2,subj][l,r]));
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
