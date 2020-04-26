function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, degF::degF{1,:H1}, m::mesh, J, dJ, phi, kubPoints, kubWeights, compoundData::compoundData{:HexToTris})
    fill!(assembledPhi,ones(Float64,1,12));
end

function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, degF::degF{1,:H1}, m::mesh, J, dJ, phi, kubPoints, kubWeights, compoundData::compoundData{:HexToKites})
    fill!(assembledPhi,ones(Float64,1,6));
end


function assemblePhi!(assembledPhi::Array{Array{Float64,2},1}, subcoord::Array{Array{Float64,2},1}, degF::degF{2,:H1div}, m::mesh, J, ddJ, jphi, kubPoints, kubWeights, compoundData::compoundData{:HexToTris})
    phi=degF.phi;
    divphi=degF.divphi;
    mt=m.meshType;
    sk=size(kubWeights);

    nPhiSubCell=size(phi,2);
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


    for subCell in 1:12
        tangent=(subcoord[subCell][:,2].-subcoord[subCell][:,1]);
        tangent=tangent./norm(tangent,2);
        jacobi!(J,ddJ,jphi,kubPoints,phi,subcoord[subCell],mt);
        for i in 1:nPhiSubCell
            for r in 1:sk[2]
                for l in 1:sk[1]
                    dotp=0.0;
                    for d in 1:m.geometry.dim
                        dotp+=jphi[d,i][l,r]*tangent[d];
                    end
                    #apply third constraint
                    A[end, (i-1)*12+subCell]+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*dotp;
                    #saving ddJ for divergence-constraint
                    #TODO: check where to compute divergence for kites since its non-constant
                    divCoeff[subCell]=copy(ddJ[1]);
                end
            end
        end
    end
    #apply second constraint (divergence)
    for i in 1:11
        #FEMEuler divergences
        A[24+i,[1,12+1,24+1]]=divCoeff[1]*[1.0,1.0,-1.0];
        A[24+i,[1+i,12+1+i,24+1+i]]=divCoeff[i+1]*[-1.0,-1.0,1.0];
        #MelvinThuburn divergences (common divergences)
    #    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    #    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
    end

    betas= A\b;

    betas2times=Array{Float64,2}(undef,nPhiSubCell,24);
    betas2times[1,1:12]  = betas[1:12];
    betas2times[1,13:24] = betas[1:12];

    betas2times[2,1:12]  = betas[13:24];
    betas2times[2,13:24] = betas[13:24];

    betas2times[3,1:12]  = betas[25:36];
    betas2times[3,13:24] = betas[25:36];

    assembledPhi[1] = betas2times[:,5:16];
    assembledPhi[2] = betas2times[:,3:14];
    assembledPhi[3] = betas2times[:,13:24];
    #minus needed for correct normals
    assembledPhi[4] = -betas2times[:,11:22];
    assembledPhi[5] = -betas2times[:,9:20];
    assembledPhi[6] = -betas2times[:,7:18];
end


function assemblePhi!(CellID::Int64, degF::degF{2,:H1div}, dJ, compoundData::compoundData{:HexToKites})
    phi=degF.phi;
    divphi=degF.divphi;
    nPhiSubCell=size(phi,2);
    nSubCells=compoundData.nSubCells;
    A=zeros(nPhiSubCell*nSubCells,nPhiSubCell*nSubCells);
    b=zeros(nPhiSubCell*nSubCells);
    tangent=zeros(2);
    divCoeff=zeros(nSubCells);
    subcoord=compoundData.getSubCells[CellID];


    #first constraint
    #outer edges
    for i in 1:nSubCells
        A[i,i]=1.0;
        b[i]=0.0;
    end
    b[1]=1.0
    b[2]=1.0

    #inner edges
    #last two subCell-ansatzfunctions correspond to inner edges
    #so second to last beta from one subCell should match last beta from next subCell
    for i in 1:nSubCells
        #add second to last beta from current subCell
        A[nSubCells+i,(nPhiSubCell-2)*nSubCells+i]=1.0;
        #add last from next subCell
        A[nSubCells+i,(nPhiSubCell-1)*nSubCells+mod(i,nSubCells)+1]=-1.0;
    #    A[12+i,24+mod(i,12)+1]=1.0;
        b[nSubCells+i]=0.0;
    end


    for subCell in 1:nSubCells
        tangent=(subcoord[subCell][:,2].-subcoord[subCell][:,1]);
        tangent=tangent./norm(tangent,2);
        jacobi!(J,ddJ,jphi,kubPoints,phi,subcoord[subCell],mt);
        for i in 1:nPhiSubCell
            for r in 1:sk[2]
                for l in 1:sk[1]
                    dotp=0.0;
                    for d in 1:m.geometry.dim
                        dotp+=jphi[d,i][l,r]*tangent[d];
                    end
                    #applying third constraint
                    A[end, (i-1)*nSubCells+subCell]+=kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*dotp;
                    #saving ddJ for divergence-constraint
                    #TODO: check where to compute divergence for kites since its non-constant
                    divCoeff[subCell]=ddJ[1];
                end
            end
        end
    end

    #apply second constraint (divergence)
    for i in 1:11
        #FEMEuler divergences
        A[+i,[1,12+1,24+1]]=[1.0,1.0,-1.0];
        A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,1.0];
        #MelvinThuburn divergences (common divergences)
    #    A[24+i,[1,12+1,24+1]]=[1.0,1.0,1.0];
    #    A[24+i,[1+i,12+1+i,24+1+i]]=[-1.0,-1.0,-1.0];
        b[24+i]=0.0;
    end

end
