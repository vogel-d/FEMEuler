function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64,
                            compoundData::compoundData)

    sk=length(quadWeights)
    kubPoints=:nothing;
    kubWeights=:nothing;

    nSubCells=compoundData.nSubCells;
    w1=zeros(sk);
    w2=zeros(sk);
    mt=m.meshType;

    center=Array{Float64,1}(undef,2);
    nSubPhiT=length(phiT);
    nSubPhiF=size(phiF,2);
    nSubPhiW=length(phiW);
    assembledPhiT1=compoundData.assembledPhi[degFT.femType];
    assembledPhiF1=compoundData.assembledPhi[degFF.femType];
    assembledPhiW1=compoundData.assembledPhi[degFW.femType];
    assembledPhiT2=compoundData.assembledPhi[degFT.femType];
    assembledPhiF2=compoundData.assembledPhi[degFF.femType];
    assembledPhiW2=compoundData.assembledPhi[degFW.femType];
    #assembledPhiT1=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF1=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW1=compoundData.assembledPhiSafe[degFW.femType];
    #assembledPhiT2=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF2=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW2=compoundData.assembledPhiSafe[degFW.femType];
    nCompoundPhiT=length(assembledPhiT1);
    nCompoundPhiF=length(assembledPhiF1);
    nCompoundPhiW=length(assembledPhiW1);

    #boundary[i] is a dict describing which subelements have an edge
    #that takes part on edge i of compound element
    compoundBoundary=compoundData.boundary;
    #assuming same number of subelements adjacent to every compound edge
    nEdgeParts=length(keys(compoundBoundary[1]));
    adjacentSubCells=Array{Int64,2}(undef,2,nEdgeParts);
    nVertices_CompoundElement=diff(m.topology.offset["20"][1:2])[1];
    coord1=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    coord2=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    subcoord1=Array{Array{Float64,2},1}(undef,nSubCells);
    subcoord2=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord1[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
        subcoord2[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    #edge integrating for compound RT0
    divphiF=degFF.divphi;
    nquadPhiF=compoundData.nquadPhi[degFF.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiF_edge=initJacobi((m.geometry.dim,size(phiF,2)),sq);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,1}(undef,sk);

    lM11=zeros(nCompoundPhiT,nCompoundPhiF);
    lM12=zeros(nCompoundPhiT,nCompoundPhiF);
    lM21=zeros(nCompoundPhiT,nCompoundPhiF);
    lM22=zeros(nCompoundPhiT,nCompoundPhiF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        adjacentSubCells1=collect(keys(compoundBoundary[eT1]));
        adjacentSubCells2=collect(keys(compoundBoundary[eT2]));
        coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord2= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, center, compoundData);
        getSubCells!(subcoord2, coord2, center, compoundData);

        assemblePhi!(assembledPhiT1, subcoord1, degFT, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF1, subcoord1, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW1, subcoord1, degFW, m, J, dJ, phiW, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiT2, subcoord2, degFT, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF2, subcoord2, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW2, subcoord2, degFW, m, J, dJ, phiW, kubPoints, kubWeights, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        if compoundData.isEdgePeriodic[e] #periodic boundary edge
            #make coordinates of current edge equal in both elements
            #has to be reversed to save correct edge direction
            #fill coord1 with Inf to avoid other coordinates being equal than the adjacent
            fill!(coord1,Inf);
            for d in 1:m.geometry.dim
                #assuming local edge i connecting local vertices i and its cyclic next neighbour in coord
                coord1[d,eT1]=coord2[d,mod(eT2,nVertices_CompoundElement)+1];
                coord1[d,mod(eT1,nVertices_CompoundElement)+1]=coord2[d,eT2];
            end
            getSubCells!(subcoord1, coord1, center, compoundData);
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
            coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
            getSubCells!(subcoord1, coord1, center, compoundData);
        else
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
        end

        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for edge_part in 1:nEdgeParts
            subCell1=adjacentSubCells[1,edge_part];
            subCell2=adjacentSubCells[2,edge_part];

            subeT1=compoundBoundary[eT1][subCell1];
            subeT2=compoundBoundary[eT2][subCell2];

            phiFn1=@views phiFtrans[subeT1];
            phiTn1=@views phiTtrans[subeT1];
            phiWn1=@views phiWtrans[subeT1];
            kubPn1=@views nquadPoints[subeT1];
            phiFn2=@views phiFtrans[subeT2];
            phiTn2=@views phiTtrans[subeT2];
            phiWn2=@views phiWtrans[subeT2];
            kubPn2=@views nquadPoints[subeT2];

            fill!(w1,0.0);
            fill!(w2,0.0);
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1+=assembledPhiW1[i][subi,subCell1]*wval[globalNumW1[i]]*phiWn1[subi];
                    @. w2+=assembledPhiW2[i][subi,subCell2]*wval[globalNumW2[i]]*phiWn2[subi];
                end
            end

            n1=m.normals[:,subeT1];
            n2=m.normals[:,subeT2];

            correctNormalsCompound!(n1,n2,eT1,eT2,compoundData);

            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    for subj in 1:nSubPhiF
                        for subi in 1:nSubPhiT
                            for r in 1:sk
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn1[subi][r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r]);
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn1[subi][r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r]);
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn2[subi][r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r]);
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn2[subi][r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r]);
                            end
                        end
                    end
                end
            end
        end

        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];

                M[gi1]+=(+0.5-gammaLoc)*lM11[i,j]*fval[gj1];
                M[gi1]+=(-0.5+gammaLoc)*lM12[i,j]*fval[gj2];
                M[gi2]+=(+0.5+gammaLoc)*lM21[i,j]*fval[gj1];
                M[gi2]+=(-0.5-gammaLoc)*lM22[i,j]*fval[gj2];
            end
        end
    end
    return nothing;
end

function discGalerkinEdges!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::Array{Float64,1}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64,
                            compoundData::compoundData)

    println("ADVSTIFFMATRIX USED")
    sk=length(quadWeights)
    kubPoints=:nothing;
    kubWeights=:nothing;

    nSubCells=compoundData.nSubCells;
    w1=zeros(sk);
    w2=zeros(sk);
    mt=m.meshType;

    center=Array{Float64,1}(undef,2);
    nSubPhiT=length(phiT);
    nSubPhiF=size(phiF,2);
    nSubPhiW=length(phiW);
    assembledPhiT1=compoundData.assembledPhi[degFT.femType];
    assembledPhiF1=compoundData.assembledPhi[degFF.femType];
    assembledPhiW1=compoundData.assembledPhi[degFW.femType];
    assembledPhiT2=compoundData.assembledPhi[degFT.femType];
    assembledPhiF2=compoundData.assembledPhi[degFF.femType];
    assembledPhiW2=compoundData.assembledPhi[degFW.femType];
    #assembledPhiT1=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF1=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW1=compoundData.assembledPhiSafe[degFW.femType];
    #assembledPhiT2=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF2=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW2=compoundData.assembledPhiSafe[degFW.femType];
    nCompoundPhiT=length(assembledPhiT1);
    nCompoundPhiF=length(assembledPhiF1);
    nCompoundPhiW=length(assembledPhiW1);

    #boundary[i] is dict describing which subelements have an edge
    #that takes part on edge i of compound element
    compoundBoundary=compoundData.boundary;
    #assuming same number of subelements adjacent to every compound edge
    nEdgeParts=length(keys(compoundBoundary[1]));
    adjacentSubCells=Array{Int64,2}(undef,2,nEdgeParts);
    nVertices_CompoundElement=diff(m.topology.offset["20"][1:2])[1];
    coord1=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    coord2=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    subcoord1=Array{Array{Float64,2},1}(undef,nSubCells);
    subcoord2=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord1[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
        subcoord2[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    #edge integrating for compound RT0
    divphiF=degFF.divphi;
    nquadPhiF=compoundData.nquadPhi[degFF.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiF_edge=initJacobi((m.geometry.dim,size(phiF,2)),sq);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,1}(undef,sk);

    lM11=zeros(nCompoundPhiT,nCompoundPhiF);
    lM12=zeros(nCompoundPhiT,nCompoundPhiF);
    lM21=zeros(nCompoundPhiT,nCompoundPhiF);
    lM22=zeros(nCompoundPhiT,nCompoundPhiF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        adjacentSubCells1=collect(keys(compoundBoundary[eT1]));
        adjacentSubCells2=collect(keys(compoundBoundary[eT2]));
        coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord2= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, center, compoundData);
        getSubCells!(subcoord2, coord2, center, compoundData);

        assemblePhi!(assembledPhiT1, subcoord1, degFT, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF1, subcoord1, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW1, subcoord1, degFW, m, J, dJ, phiW, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiT2, subcoord2, degFT, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF2, subcoord2, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW2, subcoord2, degFW, m, J, dJ, phiW, kubPoints, kubWeights, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        if compoundData.isEdgePeriodic[e] #periodic boundary edge
            #make coordinates of current edge equal in both elements
            #has to be reversed to save correct edge direction
            #fill coord1 with Inf to avoid other coordinates being equal than the adjacent
            fill!(coord1,Inf);
            for d in 1:m.geometry.dim
                #assuming local edge i connecting local vertices i and its cyclic next neighbour in coord
                coord1[d,eT1]=coord2[d,mod(eT2,nVertices_CompoundElement)+1];
                coord1[d,mod(eT1,nVertices_CompoundElement)+1]=coord2[d,eT2];
            end
            getSubCells!(subcoord1, coord1, center, compoundData);
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
            coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
            getSubCells!(subcoord1, coord1, center, compoundData);
        else
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
        end

        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for edge_part in 1:nEdgeParts
            subCell1=adjacentSubCells[1,edge_part];
            subCell2=adjacentSubCells[2,edge_part];

            subeT1=compoundBoundary[eT1][subCell1];
            subeT2=compoundBoundary[eT2][subCell2];

            phiFn1=@views phiFtrans[subeT1];
            phiTn1=@views phiTtrans[subeT1];
            phiWn1=@views phiWtrans[subeT1];
            kubPn1=@views nquadPoints[subeT1];
            phiFn2=@views phiFtrans[subeT2];
            phiTn2=@views phiTtrans[subeT2];
            phiWn2=@views phiWtrans[subeT2];
            kubPn2=@views nquadPoints[subeT2];

            fill!(w1,0.0);
            fill!(w2,0.0);
            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1+=assembledPhiW1[i][subi,subCell1]*wval[globalNumW1[i]]*phiWn1[subi];
                    @. w2+=assembledPhiW2[i][subi,subCell2]*wval[globalNumW2[i]]*phiWn2[subi];
                end
            end

            #n1=@views m.normals[:,subeT1];
            #n2=@views m.normals[:,subeT2];
            n1=m.normals[:,subeT1];
            n2=m.normals[:,subeT2];

            correctNormalsCompound!(n1,n2,eT1,eT2,compoundData);

            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    for subj in 1:nSubPhiF
                        for subi in 1:nSubPhiT
                            for r in 1:sk
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn1[subi][r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r]);
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn1[subi][r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r]);
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn2[subi][r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r]);
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn2[subi][r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r]);
                            end
                        end
                    end
                end
            end
        end

        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];

                push!(rows,gi1);
                push!(cols,gj1);
                push!(vals,(+0.5-gammaLoc)*lM11[i,j]);

                push!(rows,gi1);
                push!(cols,gj2);
                push!(vals,(-0.5+gammaLoc)*lM12[i,j]);

                push!(rows,gi2);
                push!(cols,gj1);
                push!(vals,(+0.5+gammaLoc)*lM21[i,j]);

                push!(rows,gi2);
                push!(cols,gj2);
                push!(vals,(-0.5-gammaLoc)*lM22[i,j]);
            end
        end

    end
    return nothing;
end

function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            #degFW::degF{2,S} where S,phiW::Array{Array{Float64,2},2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=length(quadWeights)

    nSubCells=compoundData.nSubCells;
    w1=[zeros(sk) for d in 1:m.geometry.dim]
    w2=[zeros(sk) for d in 1:m.geometry.dim]

    mt=m.meshType;

    center=Array{Float64,1}(undef,2);
    nSubPhiT=size(phiT,2);
    nSubPhiF=size(phiF,2);
    nSubPhiW=size(phiW,2);
    assembledPhiT1=compoundData.assembledPhi[degFT.femType];
    assembledPhiF1=compoundData.assembledPhi[degFF.femType];
    assembledPhiW1=compoundData.assembledPhi[degFW.femType];
    assembledPhiT2=compoundData.assembledPhi[degFT.femType];
    assembledPhiF2=compoundData.assembledPhi[degFF.femType];
    assembledPhiW2=compoundData.assembledPhi[degFW.femType];
    #assembledPhiT1=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF1=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW1=compoundData.assembledPhiSafe[degFW.femType];
    #assembledPhiT2=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF2=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW2=compoundData.assembledPhiSafe[degFW.femType];
    nCompoundPhiT=length(assembledPhiT1);
    nCompoundPhiF=length(assembledPhiF1);
    nCompoundPhiW=length(assembledPhiW1);

    #boundary[i] is dict describing which subelements have an edge
    #that takes part on edge i of compound element
    compoundBoundary=compoundData.boundary;
    #assuming same number of subelements adjacent to every compound edge
    nEdgeParts=length(keys(compoundBoundary[1]));
    adjacentSubCells=Array{Int64,2}(undef,2,nEdgeParts);
    nVertices_CompoundElement=diff(m.topology.offset["20"][1:2])[1];
    coord1=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    coord2=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    subcoord1=Array{Array{Float64,2},1}(undef,nSubCells);
    subcoord2=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord1[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
        subcoord2[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    J1=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ1=Array{Float64,1}(undef,sk);
    jphiWn1=initJacobi((m.geometry.dim,size(phiW,2)),sk)
    jphiTn1=initJacobi((m.geometry.dim,size(phiT,2)),sk)

    J2=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ2=Array{Float64,1}(undef,sk);
    jphiWn2=initJacobi((m.geometry.dim,size(phiW,2)),sk)
    jphiTn2=initJacobi((m.geometry.dim,size(phiT,2)),sk);

    #edge integrating for compound RT0
    divphiT=degFT.divphi;
    divphiF=degFF.divphi;
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

    lM11=zeros(nCompoundPhiT,nCompoundPhiF);
    lM12=zeros(nCompoundPhiT,nCompoundPhiF);
    lM21=zeros(nCompoundPhiT,nCompoundPhiF);
    lM22=zeros(nCompoundPhiT,nCompoundPhiF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;

        adjacentSubCells1=collect(keys(compoundBoundary[eT1]));
        adjacentSubCells2=collect(keys(compoundBoundary[eT2]));
        coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord2= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, center, compoundData);
        getSubCells!(subcoord2, coord2, center, compoundData);

        assemblePhi!(assembledPhiT1, subcoord1, m, divphiT, J_edge, ddJ_edge, jphiT_edge, nquadPhiT, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiF1, subcoord1, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW1, subcoord1, m, divphiW, J_edge, ddJ_edge, jphiW_edge, nquadPhiW, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiT2, subcoord2, m, divphiT, J_edge, ddJ_edge, jphiT_edge, nquadPhiT, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiF2, subcoord2, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW2, subcoord2, m, divphiW, J_edge, ddJ_edge, jphiW_edge, nquadPhiW, nquadPoints, quadWeights, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        if compoundData.isEdgePeriodic[e] #periodic boundary edge
            #make coordinates of current edge equal in both elements
            #has to be reversed to save correct edge direction
            #fill coord1 with Inf to avoid other coordinates being equal than the adjacent
            fill!(coord1,Inf);
            for d in 1:m.geometry.dim
                #assuming local edge i connecting local vertices i and its cyclic next neighbour in coord
                coord1[d,eT1]=                                  coord2[d,mod(eT2,nVertices_CompoundElement)+1];
                coord1[d,mod(eT1,nVertices_CompoundElement)+1]= coord2[d,eT2];
            end
            getSubCells!(subcoord1, coord1, center, compoundData);
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
            coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
            getSubCells!(subcoord1, coord1, center, compoundData);
        else
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
        end

        le=m.edgeLength[edgeData[1][e]];
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for edge_part in 1:nEdgeParts
            subCell1=adjacentSubCells[1,edge_part];
            subCell2=adjacentSubCells[2,edge_part];

            subeT1=compoundBoundary[eT1][subCell1];
            subeT2=compoundBoundary[eT2][subCell2];

            phiFn1=@views phiFtrans[subeT1];
            phiTn1=@views phiTtrans[subeT1];
            phiWn1=@views phiWtrans[subeT1];
            kubPn1=@views nquadPoints[subeT1];
            phiFn2=@views phiFtrans[subeT2];
            phiTn2=@views phiTtrans[subeT2];
            phiWn2=@views phiWtrans[subeT2];
            kubPn2=@views nquadPoints[subeT2];

            jacobi!(J1,ddJ1,jphiWn1,jphiTn1,kubPn1, phiWn1, phiTn1, subcoord1[subCell1], mt);
            jacobi!(J2,ddJ2,jphiWn2,jphiTn2,kubPn2, phiWn2, phiTn2, subcoord2[subCell2], mt);

            n1=m.normals[:,subeT1];
            n2=m.normals[:,subeT2];

            correctNormalsCompound!(n1,n2,eT1,eT2,compoundData);

            for d in 1:m.geometry.dim
                fill!(w1[d],0.0);
                fill!(w2[d],0.0);
                for i in 1:nCompoundPhiW
                    for subi in 1:nSubPhiW
                        @. w1[d]+=assembledPhiW1[i][subi,subCell1]*wval[globalNumW1[i]]*jphiWn1[d,subi];
                        @. w2[d]+=assembledPhiW2[i][subi,subCell2]*wval[globalNumW2[i]]*jphiWn2[d,subi];
                    end
                end
            end

            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    for subj in 1:nSubPhiF
                        for subi in 1:nSubPhiT
                            for r in 1:sk
                                w1jphiTn1=0.0; w2jphiTn1=0.0; w1jphiTn2=0.0; w2jphiTn2=0.0;
                                for d in 1:m.geometry.dim
                                    w1jphiTn1+=w1[d][r]*jphiTn1[d,subi][r];
                                    w2jphiTn1+=w2[d][r]*jphiTn1[d,subi][r];
                                    w1jphiTn2+=w1[d][r]*jphiTn2[d,subi][r];
                                    w2jphiTn2+=w2[d][r]*jphiTn2[d,subi][r];
                                end
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*ddJ1[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn1;
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                                           quadWeights[r]*ddJ2[r]*ddJ1[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn1;
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*ddJ2[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn2;
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                                           quadWeights[r]*ddJ2[r]*ddJ2[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn2;
                                # piola: nphiF mit 1/Je = 1/Kantenlänge, phiW und phiT mit 1/dJ*J
                            end
                        end
                    end
                end
            end
        end

        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:nCompoundPhiT
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:nCompoundPhiF
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];
                M[gi1]+=(+0.5-gammaLoc)*lM11[i,j]*fval[gj1];
                M[gi1]+=(-0.5+gammaLoc)*lM12[i,j]*fval[gj2];
                M[gi2]+=(+0.5+gammaLoc)*lM21[i,j]*fval[gj1];
                M[gi2]+=(-0.5-gammaLoc)*lM22[i,j]*fval[gj2];
            end
        end
    end
    return nothing;
end


function discGalerkinEdges!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            #degFW::degF{2,S} where S,phiW::Array{Array{Float64,2},2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=length(quadWeights)

    nSubCells=compoundData.nSubCells;
    w1=[zeros(sk) for d in 1:m.geometry.dim]
    w2=[zeros(sk) for d in 1:m.geometry.dim]

    mt=m.meshType;

    center=Array{Float64,1}(undef,m.geometry.dim);
    nSubPhiT=size(phiT,2);
    nSubPhiF=size(phiF,2);
    nSubPhiW=size(phiW,2);
    assembledPhiT1=compoundData.assembledPhi[degFT.femType];
    assembledPhiF1=compoundData.assembledPhi[degFF.femType];
    assembledPhiW1=compoundData.assembledPhi[degFW.femType];
    assembledPhiT2=compoundData.assembledPhi[degFT.femType];
    assembledPhiF2=compoundData.assembledPhi[degFF.femType];
    assembledPhiW2=compoundData.assembledPhi[degFW.femType];
    #assembledPhiT1=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF1=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW1=compoundData.assembledPhiSafe[degFW.femType];
    #assembledPhiT2=compoundData.assembledPhiSafe[degFT.femType];
    #assembledPhiF2=compoundData.assembledPhiSafe[degFF.femType];
    #assembledPhiW2=compoundData.assembledPhiSafe[degFW.femType];
    nCompoundPhiT=length(assembledPhiT1);
    nCompoundPhiF=length(assembledPhiF1);
    nCompoundPhiW=length(assembledPhiW1);

    #boundary[i] is dict describing which subelements have an edge
    #that takes part on edge i of compound element
    compoundBoundary=compoundData.boundary;
    #assuming same number of subelements adjacent to every compound edge
    nEdgeParts=length(keys(compoundBoundary[1]));
    adjacentSubCells=Array{Int64,2}(undef,2,nEdgeParts);
    nVertices_CompoundElement=diff(m.topology.offset["20"][1:2])[1];
    coord1=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    coord2=Array{Float64,2}(undef,m.geometry.dim,nVertices_CompoundElement);
    subcoord1=Array{Array{Float64,2},1}(undef,nSubCells);
    subcoord2=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing one single entry
        subcoord1[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
        subcoord2[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    J1=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ1=Array{Float64,1}(undef,sk);
    jphiWn1=initJacobi((m.geometry.dim,size(phiW,2)),sk)
    jphiTn1=initJacobi((m.geometry.dim,size(phiT,2)),sk)

    J2=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ2=Array{Float64,1}(undef,sk);
    jphiWn2=initJacobi((m.geometry.dim,size(phiW,2)),sk)
    jphiTn2=initJacobi((m.geometry.dim,size(phiT,2)),sk);

    #edge integrating for compound RT0
    divphiT=degFT.divphi;
    divphiF=degFF.divphi;
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

    lM11=zeros(nCompoundPhiT,nCompoundPhiF);
    lM12=zeros(nCompoundPhiT,nCompoundPhiF);
    lM21=zeros(nCompoundPhiT,nCompoundPhiF);
    lM22=zeros(nCompoundPhiT,nCompoundPhiF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;

        adjacentSubCells1=collect(keys(compoundBoundary[eT1]));
        adjacentSubCells2=collect(keys(compoundBoundary[eT2]));
        coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord2= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, center, compoundData);
        getSubCells!(subcoord2, coord2, center, compoundData);

        assemblePhi!(assembledPhiT1, subcoord1, m, divphiT, J_edge, ddJ_edge, jphiT_edge, nquadPhiT, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiF1, subcoord1, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW1, subcoord1, m, divphiW, J_edge, ddJ_edge, jphiW_edge, nquadPhiW, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiT2, subcoord2, m, divphiT, J_edge, ddJ_edge, jphiT_edge, nquadPhiT, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiF2, subcoord2, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        assemblePhi!(assembledPhiW2, subcoord2, m, divphiW, J_edge, ddJ_edge, jphiW_edge, nquadPhiW, nquadPoints, quadWeights, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        if compoundData.isEdgePeriodic[e] #periodic boundary edge
            #make coordinates of current edge equal in both elements
            #has to be reversed to save correct edge direction
            #fill coord1 with Inf to avoid other coordinates being equal than the adjacent
            fill!(coord1,Inf);
            for d in 1:m.geometry.dim
                #assuming local edge i connecting local vertices i and its cyclic next neighbour in coord
                coord1[d,eT1]=                                  coord2[d,mod(eT2,nVertices_CompoundElement)+1];
                coord1[d,mod(eT1,nVertices_CompoundElement)+1]= coord2[d,eT2];
            end
            getSubCells!(subcoord1, coord1, center, compoundData);
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
            coord1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
            getSubCells!(subcoord1, coord1, center, compoundData);
        else
            orderAdjacentSubCells!(adjacentSubCells,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);
        end

        le=m.edgeLength[edgeData[1][e]];
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);
        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumT2,degFT,inc2);

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for edge_part in 1:nEdgeParts
            subCell1=adjacentSubCells[1,edge_part];
            subCell2=adjacentSubCells[2,edge_part];

            subeT1=compoundBoundary[eT1][subCell1];
            subeT2=compoundBoundary[eT2][subCell2];

            phiFn1=@views phiFtrans[subeT1];
            phiTn1=@views phiTtrans[subeT1];
            phiWn1=@views phiWtrans[subeT1];
            kubPn1=@views nquadPoints[subeT1];
            phiFn2=@views phiFtrans[subeT2];
            phiTn2=@views phiTtrans[subeT2];
            phiWn2=@views phiWtrans[subeT2];
            kubPn2=@views nquadPoints[subeT2];

            jacobi!(J1,ddJ1,jphiWn1,jphiTn1,kubPn1, phiWn1, phiTn1, subcoord1[subCell1], mt);
            jacobi!(J2,ddJ2,jphiWn2,jphiTn2,kubPn2, phiWn2, phiTn2, subcoord2[subCell2], mt);

            n1=m.normals[:,subeT1];
            n2=m.normals[:,subeT2];

            correctNormalsCompound!(n1,n2,eT1,eT2,compoundData);

            for d in 1:m.geometry.dim
                fill!(w1[d],0.0);
                fill!(w2[d],0.0);
                for i in 1:nCompoundPhiW
                    for subi in 1:nSubPhiW
                        @. w1[d]+=assembledPhiW1[i][subi,subCell1]*wval[globalNumW1[i]]*jphiWn1[d,subi];
                        @. w2[d]+=assembledPhiW2[i][subi,subCell2]*wval[globalNumW2[i]]*jphiWn2[d,subi];
                    end
                end
            end

            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    for subj in 1:nSubPhiF
                        for subi in 1:nSubPhiT
                            for r in 1:sk
                                w1jphiTn1=0.0; w2jphiTn1=0.0; w1jphiTn2=0.0; w2jphiTn2=0.0;
                                for d in 1:m.geometry.dim
                                    w1jphiTn1+=w1[d][r]*jphiTn1[d,subi][r];
                                    w2jphiTn1+=w2[d][r]*jphiTn1[d,subi][r];
                                    w1jphiTn2+=w1[d][r]*jphiTn2[d,subi][r];
                                    w2jphiTn2+=w2[d][r]*jphiTn2[d,subi][r];
                                end
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*ddJ1[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn1;
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                                           quadWeights[r]*ddJ2[r]*ddJ1[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn1;
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*ddJ2[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn2;
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                                           quadWeights[r]*ddJ2[r]*ddJ2[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn2;
                                # piola: nphiF mit 1/Je = 1/Kantenlänge, phiW und phiT mit 1/dJ*J
                            end
                        end
                    end
                end
            end
        end

        for i in 1:nCompoundPhiT
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:nCompoundPhiF
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];

                push!(rows,gi1);
                push!(cols,gj1);
                push!(vals,(+0.5-gammaLoc)*lM11[i,j]);

                push!(rows,gi1);
                push!(cols,gj2);
                push!(vals,(-0.5+gammaLoc)*lM12[i,j]);

                push!(rows,gi2);
                push!(cols,gj1);
                push!(vals,(+0.5+gammaLoc)*lM21[i,j]);

                push!(rows,gi2);
                push!(cols,gj2);
                push!(vals,(-0.5-gammaLoc)*lM22[i,j]);
            end
        end
    end
    return nothing;
end



function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{2,:H1xH1},phiW::Array{Array{Float64,2},2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64, coord::Array{Float64,2},
                            compoundData::compoundData)

    sk=length(quadWeights)

    nSubCells=compoundData.nSubCells;
    w_clean=zeros(sk);
    w1=[zeros(sk) for d in 1:m.geometry.dim]
    w2=[zeros(sk) for d in 1:m.geometry.dim]

    mt=m.meshType;

    nSubPhiW=length(phiW);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    #boundary[i] is dict describing which subelements have an edge
    #that takes part on edge i of compound element
    compoundBoundary=compoundData.boundary;
    #assuming same number of subelements adjacent to every compound edge
    nEdgeParts=length(keys(compoundBoundary[1]));

    J1=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ1=Array{Float64,1}(undef,sk);
    jphiTn1=initJacobi((m.geometry.dim,size(phiT,2)),sk)

    J2=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ2=Array{Float64,1}(undef,sk);
    jphiTn2=initJacobi((m.geometry.dim,size(phiT,2)),sk);

    w1=[zeros(sk) for d in 1:m.geometry.dim]
    w2=[zeros(sk) for d in 1:m.geometry.dim]

    lM11=zeros(nT,nF);
    lM12=zeros(nT,nF);
    lM21=zeros(nT,nF);
    lM22=zeros(nT,nF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;

        adjacentSubCells1=collect(keys(compoundBoundary[eT1]));
        adjacentSubCells2=collect(keys(compoundBoundary[eT2]));
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, compoundData);
        getSubCells!(subcoord2, coord2, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        orderAdjacentSubCells!(adjacentSubCells,m,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);

        le=m.edgeLength[edgeData[1][e]];
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);

        for edge_part in 1:nEdgeParts
            subCell1=adjacentSubCells[1,edge_part];
            subCell2=adjacentSubCells[2,edge_part];

            jacobi!(J1,ddJ1,jphiTn1,kubPn1,phiTn1,subcoord1[subCell1],mt);
            jacobi!(J2,ddJ2,jphiTn2,kubPn2,phiTn2,subcoord1[subCell1],mt);

            subeT1=compoundBoundary[eT1][subCell1];
            subeT2=compoundBoundary[eT2][subCell2];

            n1=@views m.normals[:,subeT1];
            n2=@views m.normals[:,subeT2];
            le=(1/nEdgeParts)*m.edgeLength[edgeData[1][e]];


            for d in 1:m.geometry.dim
                fill!(w1[d],w_clean);
                fill!(w2[d],w_clean);
                for i in 1:nCompoundPhiW
                    for subi in 1:nSubPhiW
                        @. w1[d]+=assembledPhi[i][subi,subCell1]*wval[globalNumW1[subi]]*phiWn1[subi];
                        @. w2[d]+=assembledPhi[i][subi,subCell2]*wval[globalNumW2[subi]]*phiWn2[subi];
                    end
                end
            end


            phiFn1=@views phiFtrans[subeT1];
            phiTn1=@views phiTtrans[subeT1];
            phiWn1=@views phiWtrans[subeT1];
            kubPn1=@views nquadPoints[subeT1];
            phiFn2=@views phiFtrans[subeT2];
            phiTn2=@views phiTtrans[subeT2];
            phiWn2=@views phiWtrans[subeT2];
            kubPn2=@views nquadPoints[subeT2];

            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    for subj in nSubPhiF
                        for subi in nSubPhiT
                            for r in 1:sk
                                w1jphiTn1=0.0; w2jphiTn1=0.0; w1jphiTn2=0.0; w2jphiTn2=0.0;
                                for d in 1:m.geometry.dim
                                    w1jphiTn1+=w1[d][r]*jphiTn1[d,subi][r];
                                    w2jphiTn1+=w2[d][r]*jphiTn1[d,subi][r];
                                    w1jphiTn2+=w1[d][r]*jphiTn2[d,subi][r];
                                    w2jphiTn2+=w2[d][r]*jphiTn2[d,subi][r];
                                end
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn1;
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn1;
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ2[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn2;
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell1]*assembledPhiF2[j][subj,subCell1]*
                                           quadWeights[r]*ddJ2[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn2;
                                           # piola: nphiF mit 1/Je = 1/Kantenlänge, w nicht transformiert wg H1xH1, phiT mit 1/dJ*J
                            end
                        end
                    end
                end
            end
        end

        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];
                M[gi1]+=(+0.5-gammaLoc)*lM11[i,j]*fval[gj1];
                M[gi1]+=(-0.5+gammaLoc)*lM12[i,j]*fval[gj2];
                M[gi2]+=(+0.5+gammaLoc)*lM21[i,j]*fval[gj1];
                M[gi2]+=(-0.5-gammaLoc)*lM22[i,j]*fval[gj2];
            end
        end
    end
    return nothing;
end















#=
for j in 1:nCompoundPhiF
    for i in 1:nCompoundPhiT
        for subj in 1:nSubPhiF
            for subi in 1:nSubPhiT
                for r in 1:sk
                    w1jphiTn1=0.0; w2jphiTn1=0.0; w1jphiTn2=0.0; w2jphiTn2=0.0;
                    for d in 1:m.geometry.dim
                        w1jphiTn1+=w1[d][r]*jphiTn1[d,subi][r];
                        w2jphiTn1+=w2[d][r]*jphiTn1[d,subi][r];
                        w1jphiTn2+=w1[d][r]*jphiTn2[d,subi][r];
                        w2jphiTn2+=w2[d][r]*jphiTn2[d,subi][r];
                    end
                    lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                               quadWeights[r]*ddJ1[r]*ddJ1[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn1;
                    lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                               quadWeights[r]*ddJ2[r]*ddJ1[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn1;
                    lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                               quadWeights[r]*ddJ1[r]*ddJ2[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn2;
                    lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                               quadWeights[r]*ddJ2[r]*ddJ2[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn2;
                    # piola: nphiF mit 1/Je = 1/Kantenlänge, phiW und phiT mit 1/dJ*J
                end
            end
        end
        if e==3 && edge_part==1
            if i==1 && (j==2 || j==4)
                #println(w1jphiTn1)
                #println(w2jphiTn1)
                println(i)
                println(j)
                println(lM11[i,j])
                println(lM12[i,j])
            end
        end
    end
end
=#
