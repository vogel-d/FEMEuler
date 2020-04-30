function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64,
                            compoundData::compoundData)

    sk=length(quadWeights)

    nSubCells=compoundData.nSubCells;
    w_clean=zeros(sk);
    w1=Array{Array{Float64,2},1}(undef,nSubCells);
    w2=Array{Array{Float64,2},1}(undef,nSubCells);
    mt=m.meshType;

    nSubPhiW=length(phiW);
    assembledPhiT=compoundData.assembledPhi[degFT.femType];
    assembledPhiF=compoundData.assembledPhi[degFF.femType];
    assembledPhiW=compoundData.assembledPhi[degFW.femType];
    nCompoundPhiT=length(assembledPhiT);
    nCompoundPhiF=length(assembledPhiF);
    nCompoundPhiW=length(assembledPhiW);

    #boundary[i] is a dict describing which subelements have an edge
    #that takes part on edge i of compound element
    compoundBoundary=compoundData.boundary;
    #assuming same number of subelements adjacent to every compound edge
    nEdgeParts=length(keys(compoundBoundary[1]));
    adjacentSubCells=Array{Int64,2}(undef,2,nEdgeParts);
    coord1=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    coord2=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    subcoord1=Array{Array{Float64,2},1}(undef,nSubCells);
    subcoord2=Array{Array{Float64,2},1}(undef,nSubCells);

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
        adjacentSubCells2=collect(keys(compoundBboundary[eT2]));
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, compoundData);
        getSubCells!(subcoord2, coord2, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        orderAdjacentSubCells!(adjacentSubCells,m,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);

        fill!(w1,w_clean);
        fill!(w2,w_clean);
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

            subeT1=boundary[eT1][subCell1];
            subeT2=boundary[eT2][subCell2];

            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1+=assembledPhi[i][subi,subCell1]*wval[globalNumW1[subi]]*phiWn1[subi];
                    @. w2+=assembledPhi[i][subi,subCell2]*wval[globalNumW2[subi]]*phiWn2[subi];
                end
            end

            n1=@views m.normals[:,subeT1];
            n2=@views m.normals[:,subeT2];

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
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn1[i][r]*(n1[1]*phiFn1[1,j][r]+n1[2]*phiFn1[2,j][r]);
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn1[i][r]*(n2[1]*phiFn2[1,j][r]+n2[2]*phiFn2[2,j][r]);
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn2[i][r]*(n1[1]*phiFn1[1,j][r]+n1[2]*phiFn1[2,j][r]);
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn2[i][r]*(n2[1]*phiFn2[1,j][r]+n2[2]*phiFn2[2,j][r]);
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


    sk=length(quadWeights)

    nSubCells=compoundData.nSubCells;
    w_clean=zeros(sk);
    w1=zeros(sk);
    w2=zeros(sk);
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
    adjacentSubCells=Array{Int64,2}(undef,2,nEdgeParts);
    coord1=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    coord2=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    subcoord1=Array{Array{Float64,2},1}(undef,nSubCells);
    subcoord2=Array{Array{Float64,2},1}(undef,nSubCells);


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
        adjacentSubCells2=collect(keys(compoundBboundary[eT2]));
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
        getSubCells!(subcoord1, coord1, compoundData);
        getSubCells!(subcoord2, coord2, compoundData);

        #order adjacentSubCells to have adjacentSubCells[1,j] sharing an edge with adjacentSubCells[2,j]
        orderAdjacentSubCells!(adjacentSubCells,m,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2);

        fill!(w1,w_clean);
        fill!(w2,w_clean);
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

            subeT1=boundary[eT1][subCell1];
            subeT2=boundary[eT2][subCell2];

            for i in 1:nCompoundPhiW
                for subi in 1:nSubPhiW
                    @. w1+=assembledPhi[i][subi,subCell1]*wval[globalNumW1[subi]]*phiWn1[subi];
                    @. w2+=assembledPhi[i][subi,subCell2]*wval[globalNumW2[subi]]*phiWn2[subi];
                end
            end

            n1=@views m.normals[:,subeT1];
            n2=@views m.normals[:,subeT2];

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
                                lM11[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn1[i][r]*(n1[1]*phiFn1[1,j][r]+n1[2]*phiFn1[2,j][r]);
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn1[i][r]*(n2[1]*phiFn2[1,j][r]+n2[2]*phiFn2[2,j][r]);
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF1[j][subj,subCell1]*
                                            quadWeights[r]*w1[r]*phiTn2[i][r]*(n1[1]*phiFn1[1,j][r]+n1[2]*phiFn1[2,j][r]);
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell2]*assembledPhiF2[j][subj,subCell2]*
                                            quadWeights[r]*w2[r]*phiTn2[i][r]*(n2[1]*phiFn2[1,j][r]+n2[2]*phiFn2[2,j][r]);
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
    jphiWn1=initJacobi((m.geometry.dim,size(phiW,2)),sk)
    jphiTn1=initJacobi((m.geometry.dim,size(phiT,2)),sk)

    J2=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ2=Array{Float64,1}(undef,sk);
    jphiWn2=initJacobi((m.geometry.dim,size(phiW,2)),sk)
    jphiTn2=initJacobi((m.geometry.dim,size(phiT,2)),sk);

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
        adjacentSubCells2=collect(keys(compoundBboundary[eT2]));
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

            jacobi!(J1,ddJ1,jphiWn1,jphiTn1,kubPn1, phiWn1, phiTn1, subcoord1[subCell1], mt);
            jacobi!(J2,ddJ2,jphiWn2,jphiTn2,kubPn2, phiWn2, phiTn2, subcoord2[subCell2], mt);

            subeT1=boundary[eT1][subCell1];
            subeT2=boundary[eT2][subCell2];

            n1=@views m.normals[:,subeT1];
            n2=@views m.normals[:,subeT2];
            le=(1/nEdgeParts)*m.edgeLength[edgeData[1][e]];


            for d in 1:m.geometry.dim
                fill!(w1[d],w_clean);
                fill!(w2[d],w_clean);
                for i in 1:nCompoundPhiW
                    for subi in 1:nSubPhiW
                        @. w1[d]+=assembledPhi[i][subi,subCell1]*wval[globalNumW1[subi]]*jphiWn1[subi];
                        @. w2[d]+=assembledPhi[i][subi,subCell2]*wval[globalNumW2[subi]]*jphiWn2[subi];
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
                                           quadWeights[r]*ddJ1[r]*ddJ1[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn1;
                                lM12[i,j]+=assembledPhiT1[i][subi,subCell1]*assembledPhiF2[j][subj,subCell1]*
                                           quadWeights[r]*ddJ2[r]*ddJ1[r]*(n2[1]*phiFn2[1,subj][r]+n2[2]*phiFn2[2,subj][r])*w2jphiTn1;
                                lM21[i,j]+=assembledPhiT2[i][subi,subCell1]*assembledPhiF1[j][subj,subCell1]*
                                           quadWeights[r]*ddJ1[r]*ddJ2[r]*(n1[1]*phiFn1[1,subj][r]+n1[2]*phiFn1[2,subj][r])*w1jphiTn2;
                                lM22[i,j]+=assembledPhiT2[i][subi,subCell1]*assembledPhiF2[j][subj,subCell1]*
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
        adjacentSubCells2=collect(keys(compoundBboundary[eT2]));
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

            subeT1=boundary[eT1][subCell1];
            subeT2=boundary[eT2][subCell2];

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
