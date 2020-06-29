function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, compoundData::compoundData{:HexToTris})
    fill!(assembledPhi,ones(Float64,1,12));
end

function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, compoundData::compoundData{:HexToKites})
    fill!(assembledPhi,ones(Float64,1,6));
end

function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, compoundData::compoundData{:RectToKites})
    fill!(assembledPhi,ones(Float64,1,4));
end

function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, compoundData::compoundData{:RectToTris})
    fill!(assembledPhi,ones(Float64,1,8));
end


function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, m::mesh, divphi, J, ddJ, jphi, nquadPhi, nquadPoints, quadWeights, compoundData::compoundData{:HexToTris})
    mt=m.meshType;

    nPhiSubCell=size(jphi,2);
    nSubCells=compoundData.nSubCells;
    A=zeros(nPhiSubCell*nSubCells,nPhiSubCell*nSubCells);
    b=zeros(nPhiSubCell*nSubCells);
    tangent=zeros(2);
    divCoeff=zeros(nSubCells);

    #first constraint
    #outer edges
    for i in 1:12
        A[i,i]=1.0;
        b[i]=0.0;
    end
    b[1]=1.0
    b[2]=1.0

    #inner edges
    #last two subCell-ansatzfunctions correspond to inner edges
    #so second to last beta from one subCell should match last beta from next subCell
    for i in 1:12
        #add second to last beta from current subCell
        A[12+i,12+i]=1.0;
        #add last from next subCell
        A[12+i,24+mod(i,12)+1]=-1.0;
        #A[12+i,24+mod(i,12)+1]=1.0;
    end

    dotp=0.0;
    for subCell in 1:12
        for i in 1:2
            tangent[i]=subcoord[subCell][i,2]-subcoord[subCell][i,1];
        end
        @. tangent=tangent/sqrt(tangent[1]^2+tangent[2]^2);
        jacobi!(J,ddJ,jphi,nquadPoints[1],nquadPhi[1],subcoord[subCell],mt);
        for i in 1:nPhiSubCell
            for r in 1:length(quadWeights)
                dotp=0.0;
                for d in 1:m.geometry.dim
                    dotp+=jphi[d,i][r]*tangent[d];
                end
                #apply third constraint
                A[end, (i-1)*12+subCell]+=quadWeights[r]*(ddJ[r]/abs(ddJ[r]))*dotp;
            end
            #saving ddJ for divergence-constraint
            divCoeff[subCell]=copy(ddJ[1]);
        end
    end

    #apply second constraint (divergence)
    for i in 1:11
        #FEMEuler divergences
        A[24+i,1]=      1*divCoeff[1];
        A[24+i,12+1]=   1*divCoeff[1];
        A[24+i,24+1]=   (-1)*divCoeff[1];

        A[24+i,1+i]=    (-1)*divCoeff[i+1];
        A[24+i,12+1+i]= (-1)*divCoeff[i+1];
        A[24+i,24+1+i]= 1*divCoeff[i+1];
        #MelvinThuburn divergences (common divergences)
    #    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    #    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    end

    betas= A\b;

    for compoundPhi in 1:3
        for subPhi in 1:3
            for subElement in 1:12
                #                                                                   walking through 1:12 cyclic; starting from 1, then 11, then 9, ...
                assembledPhi[compoundPhi][subPhi,subElement] = betas[(subPhi-1)*12+(1+mod(-1+subElement-2*(compoundPhi-1),12))]
            end
        end
    end
    for compoundPhi in 4:6
        for subPhi in 1:3
            for subElement in 1:12
                #                                                                   walking through 1:12 cyclic; starting from 1, then 11, then 9, ...
                assembledPhi[compoundPhi][subPhi,subElement] = -betas[(subPhi-1)*12+(1+mod(-1+subElement-2*(compoundPhi-1),12))]
            end
        end
    end
end



function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, m::mesh, divphi, J, ddJ, jphi, nquadPhi, nquadPoints, quadWeights, compoundData::compoundData{:HexToKites})
    mt=m.meshType;

    nPhiSubCell=size(jphi,2);
    nSubCells=compoundData.nSubCells;
    A=Array{Float64,2}(undef,nPhiSubCell*nSubCells,nPhiSubCell*nSubCells);
    fill!(A,0.0)
    b=Array{Float64,1}(undef,nPhiSubCell*nSubCells);
    fill!(b,0.0)
    tangent=Array{Float64,1}(undef,2);
    divCoeff=Array{Float64,1}(undef,nSubCells);

    #first constraint
    #outer edges
    for i in 1:12
        A[i,i]=1.0;
        b[i]=0.0;
    end
    b[7]=1.0;
    b[2]=1.0;

    #inner edges
    #last two subCell-ansatzfunctions correspond to inner edges
    #so second to last beta from one subCell should match last beta from next subCell
    for i in 1:6
        #add second to last beta from current subCell
        A[12+i,12+i]=1.0;
        #add last from next subCell
        #A[12+i,18+mod(i,6)+1]=-1.0;
        A[12+i,18+mod(i,6)+1]=1.0;
    end

    dotp=0.0;
    for subCell in 1:6
        for edge in 1:2
            for i in 1:2
                tangent[i]=subcoord[subCell][i,edge+1]-subcoord[subCell][i,edge];
            end
            @. tangent=tangent/sqrt(tangent[1]^2+tangent[2]^2);
            jacobi!(J,ddJ,jphi,nquadPoints[edge],nquadPhi[edge],subcoord[subCell],mt);
            for i in 1:nPhiSubCell
                for r in 1:length(quadWeights)
                    dotp=0.0;
                    for d in 1:m.geometry.dim
                        dotp+=jphi[d,i][r]*tangent[d];
                    end
                    #apply third constraint
                    A[end, (i-1)*6+subCell]+=quadWeights[r]*(ddJ[r]/abs(ddJ[r]))*dotp;
                    #saving ddJ for divergence-constraint
                    divCoeff[subCell]=copy(ddJ[1]);
                    if divCoeff[subCell]<0.0
                        println("negative determinant while assembling phi. not considered");
                    end
                end
            end
        end
    end

    #TODO: check where to compute divergence since its non-constant
    divCoeff=ones(nSubCells);
    #apply second constraint (divergence)
    for i in 1:5
        #FEMEuler divergences
        A[18+i,1]=      1*divCoeff[1];
        A[18+i,6+1]=    1*divCoeff[1];
        A[18+i,12+1]=   (-1)*divCoeff[1];
        A[18+i,18+1]=   (-1)*divCoeff[1];

        A[18+i,1+i]=    (-1)*divCoeff[i+1];
        A[18+i,6+1+i]=  (-1)*divCoeff[i+1];
        A[18+i,12+1+i]= 1*divCoeff[i+1];
        A[18+i,18+1+i]= 1*divCoeff[i+1];
        #MelvinThuburn divergences (common divergences)
    #    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    #    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    end

    betas= A\b;

    for compoundPhi in 1:3
        for subPhi in 1:4
            for subElement in 1:6
                #                                                                   walking through 1:6 cyclic; starting from 1, then 6, then 5, ...
                assembledPhi[compoundPhi][subPhi,subElement] = betas[(subPhi-1)*6+(1+mod(-1+subElement-(compoundPhi-1),6))]
            end
        end
    end
    for compoundPhi in 4:6
        for subPhi in 1:4
            for subElement in 1:6
                #                                                                   walking through 1:6 cyclic; starting from 1, then 6, then 5, ...
                assembledPhi[compoundPhi][subPhi,subElement] = -betas[(subPhi-1)*6+(1+mod(-1+subElement-(compoundPhi-1),6))]
            end
        end
    end
end


function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, m::mesh, divphi, J, ddJ, jphi, nquadPhi, nquadPoints, quadWeights, compoundData::compoundData{:RectToKites})
    mt=m.meshType;

    nPhiSubCell=size(jphi,2);
    nSubCells=compoundData.nSubCells;
    A=Array{Float64,2}(undef,nPhiSubCell*nSubCells,nPhiSubCell*nSubCells);
    fill!(A,0.0)
    b=Array{Float64,1}(undef,nPhiSubCell*nSubCells);
    fill!(b,0.0)
    tangent=Array{Float64,1}(undef,2);
    divCoeff=Array{Float64,1}(undef,nSubCells);

    #first constraint
    #outer edges
    for i in 1:8
        A[i,i]=1.0;
    end
    b[5]=1.0;
    b[2]=1.0;

    #inner edges
    #last two subCell-ansatzfunctions correspond to inner edges
    #so second to last beta from one subCell should match last beta from next subCell
    for i in 1:4
        #add second to last beta from current subCell
        A[8+i,8+i]=1.0;
        #add last from next subCell
        A[8+i,12+mod(i,4)+1]=1.0;
    end

    dotp=0.0;
    for subCell in 1:4
        for edge in 1:2
            for i in 1:2
                tangent[i]=subcoord[subCell][i,edge+1]-subcoord[subCell][i,edge];
            end
            @. tangent=tangent/sqrt(tangent[1]^2+tangent[2]^2);
            jacobi!(J,ddJ,jphi,nquadPoints[edge],nquadPhi[edge],subcoord[subCell],mt);
            for i in 1:nPhiSubCell
                for r in 1:length(quadWeights)
                    dotp=0.0;
                    for d in 1:m.geometry.dim
                        dotp+=jphi[d,i][r]*tangent[d];
                    end
                    #apply third constraint
                    A[end, (i-1)*4+subCell]+=quadWeights[r]*(ddJ[r]/abs(ddJ[r]))*dotp;
                    #saving ddJ for divergence-constraint
                    divCoeff[subCell]=copy(ddJ[1]);
                    if divCoeff[subCell]<0.0
                        println("negative determinant while assembling phi. not considered");
                    end
                end
            end
        end
    end

    #TODO: check where to compute divergence since its non-constant
    divCoeff=ones(nSubCells);
    #apply second constraint (divergence)
    for i in 1:3
        #FEMEuler divergences
        A[12+i,1]=      1*divCoeff[1];
        A[12+i,4+1]=    1*divCoeff[1];
        A[12+i,8+1]=    (-1)*divCoeff[1];
        A[12+i,12+1]=   (-1)*divCoeff[1];

        A[12+i,1+i]=    (-1)*divCoeff[i+1];
        A[12+i,4+1+i]=  (-1)*divCoeff[i+1];
        A[12+i,8+1+i]=  1*divCoeff[i+1];
        A[12+i,12+1+i]= 1*divCoeff[i+1];
        #MelvinThuburn divergences (common divergences)
    #    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    #    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    end

    betas= A\b;

    for compoundPhi in 1:2
        for subPhi in 1:nPhiSubCell
            for subElement in 1:nSubCells
                #                                                                   walking through 1:4 cyclic; starting from 1, then 4, then 3, ...
                assembledPhi[compoundPhi][subPhi,subElement] = betas[(subPhi-1)*4+(1+mod(-1+subElement-(compoundPhi-1),4))]
            end
        end
    end
    for compoundPhi in 3:4
        for subPhi in 1:nPhiSubCell
            for subElement in 1:nSubCells
                #                                                                   walking through 1:4 cyclic; starting from 1, then 4, then 3, ...
                assembledPhi[compoundPhi][subPhi,subElement] = -betas[(subPhi-1)*4+(1+mod(-1+subElement-(compoundPhi-1),4))]
            end
        end
    end
end


function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, m::mesh, divphi, J, ddJ, jphi, nquadPhi, nquadPoints, quadWeights, compoundData::compoundData{:RectToTris})
    mt=m.meshType;

    nPhiSubCell=size(jphi,2);
    nSubCells=compoundData.nSubCells;
    A=Array{Float64,2}(undef,nPhiSubCell*nSubCells,nPhiSubCell*nSubCells);
    fill!(A,0.0)
    b=Array{Float64,1}(undef,nPhiSubCell*nSubCells);
    fill!(b,0.0)
    tangent=Array{Float64,1}(undef,2);
    divCoeff=Array{Float64,1}(undef,nSubCells);

    #first constraint
    #outer edges
    for i in 1:8
        A[i,i]=1.0;
    end
    b[1]=1.0;
    b[2]=1.0;

    #inner edges
    #last two subCell-ansatzfunctions correspond to inner edges
    #so second to last beta from one subCell should match last beta from next subCell
    for i in 1:8
        #add second to last beta from current subCell
        A[8+i,8+i]=1.0;
        #add last from next subCell
        A[8+i,16+mod(i,8)+1]=-1.0;
    end

    dotp=0.0;
    for subCell in 1:8
        for i in 1:2
            tangent[i]=subcoord[subCell][i,2]-subcoord[subCell][i,1];
        end
        @. tangent=tangent/sqrt(tangent[1]^2+tangent[2]^2);
        jacobi!(J,ddJ,jphi,nquadPoints[1],nquadPhi[1],subcoord[subCell],mt);
        for i in 1:nPhiSubCell
            for r in 1:length(quadWeights)
                dotp=0.0;
                for d in 1:m.geometry.dim
                    dotp+=jphi[d,i][r]*tangent[d];
                end
                #apply third constraint
                A[end, (i-1)*8+subCell]+=quadWeights[r]*(ddJ[r]/abs(ddJ[r]))*dotp;
                #saving ddJ for divergence-constraint
                divCoeff[subCell]=copy(ddJ[1]);
                if divCoeff[subCell]<0.0
                    println("negative determinant while assembling phi. not considered");
                end
            end
        end
    end
    #println(A[24,:])

    #TODO: check where to compute divergence since its non-constant
    divCoeff=ones(nSubCells);
    #apply second constraint (divergence)
    for i in 1:7
        #FEMEuler divergences
        A[16+i,1]=      1*divCoeff[1];
        A[16+i,8+1]=    1*divCoeff[1];
        A[16+i,16+1]=    (-1)*divCoeff[1];

        A[16+i,1+i]=    (-1)*divCoeff[i+1];
        A[16+i,8+1+i]=  (-1)*divCoeff[i+1];
        A[16+i,16+1+i]=  1*divCoeff[i+1];
        #MelvinThuburn divergences (common divergences)
    #    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    #    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    end

    betas= A\b;
global betaz=betas
global Az=A
    for compoundPhi in 1:2
        for subPhi in 1:nPhiSubCell
            for subElement in 1:nSubCells
                #                                                                   walking through 1:8 cyclic; starting from 1, then 11, then 9, ...
                assembledPhi[compoundPhi][subPhi,subElement] = betas[(subPhi-1)*8+(1+mod(-1+subElement-2*(compoundPhi-1),8))]
                #assembledPhi[compoundPhi][subPhi,subElement] = betas[(subPhi-1)*8+]
            end
        end
    end
    for compoundPhi in 3:4
        for subPhi in 1:nPhiSubCell
            for subElement in 1:nSubCells
                #                                                                   walking through 1:8 cyclic; starting from 1, then 11, then 9, ...
                assembledPhi[compoundPhi][subPhi,subElement] = -betas[(subPhi-1)*8+(1+mod(-1+subElement-2*(compoundPhi-1),8))]
            end
        end
    end
end



function assemblePhiHexToKites(femElements::Set{Symbol})
    assembledPhi=Dict();
    for type in femElements
        if type==:RT0
            #phi2times includes the factors one compound ansatzfunction needs
            #from each subelement (coloumns) with its ansatzfunctions (rows)
            #assemblation of every individual phi shows cyclic behaviour
            phi2times = (1/3)*
                        [0 3  0  0  0  0 0 3  0  0  0  0;
                         3 0  0  0  0  0 3 0  0  0  0  0;
                         0 2  1  0 -1 -2 0 2  1  0 -1 -2;
                         2 0 -2 -1  0  1 2 0 -2 -1  0  1];

            phi1 = phi2times[:,7:12];
            phi2 = phi2times[:,6:11];
            phi3 = phi2times[:,5:10];
            phi4 = -phi2times[:,4:9];
            phi5 = -phi2times[:,3:8];
            phi6 = -phi2times[:,2:7];

            assembledPhi[type]=[phi1, phi2, phi3, phi4, phi5, phi6];
        elseif type==:DG0
            phi2times = ones(Float64,1,12);

            phi1 = phi2times[:,7:12];

            assembledPhi[type]=[phi1];
        else
            error("Tell me how to assemble HexToKite compound ansatzfunctions in $type")
        end
    end
    return assembledPhi;
end
