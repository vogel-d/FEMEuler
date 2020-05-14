function assembStiffCompound!(p::femProblem)
    if p.type==:boussinesq
        degF=p.degFBoundary;
        pkey=p.femType[:p][1];
        bkey=p.femType[:b][1];
        vkey=p.femType[:v][1];

        z=[0.0, 1.0];

        Spv=assembStiffCompound(degF[pkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints, p.compoundData);
        Svp = copy(-Spv');
        Sbv=assembStiffCompound(degF[bkey], degF[vkey], z, p.mesh, p.kubWeights, p.kubPoints, p.compoundData)
        Svb = copy(Sbv');

        p.stiffM[:pv]=Spv[1:degF[pkey].num,:];
        p.stiffM[:bv]=Sbv[1:degF[bkey].num,:];
        p.stiffM[:vp]=Svp[1:degF[vkey].num,:];
        p.stiffM[:vb]=Svb[1:degF[vkey].num,:];

    elseif p.type==:compressible
        degF=p.degFBoundary;
        rhokey=p.femType[:rho][1];
        vkey=p.femType[:rhoV][1];
        pkey=p.femType[:p][1];

        z=[0.0,1.0]

        Spv=assembStiffCompound(degF[pkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints, p.compoundData);
        Svp = copy(-Spv');
        Srhov=assembStiffCompound(degF[rhokey], degF[vkey], z, p.mesh, p.kubWeights, p.kubPoints, p.compoundData)
        Svrho = copy(Srhov');

        p.stiffM[:rho]=Spv[1:degF[rhokey].num,:];
        p.stiffM[:vp]=Svp[1:degF[vkey].num,:];
        p.stiffM[:vrho]=Svrho[1:degF[vkey].num,:];

    elseif p.type==:shallow
        degF=p.degFBoundary;
        rhokey=p.femType[:rho][1];
        vkey=p.femType[:rhoV][1];
        pkey=p.femType[:p][1];

        Spv=assembStiffCompound(degF[pkey], degF[vkey], p.mesh, p.kubWeights, p.kubPoints, p.compoundData);
        Svp = copy(-Spv');

        p.stiffM[:rho]=Spv[1:degF[rhokey].num,:];
        p.stiffM[:vp]=Svp[1:degF[vkey].num,:];

    else
        comp=Set{Symbol}()
        for i in collect(keys(p.femType))
            push!(comp,p.femType[i][1])
        end
        for i in comp
            assembStiffCompound!(p,i);
        end
    end
end

#skalare Größe mit Divergenz von vektorieller Größe
function assembStiffCompound(degFs::degF{1,:H1}, degFv::degF{2,:H1div}, m::mesh, kubWeights::Array{Float64,2}, kubPoints::Array{Float64,2}, compoundData::compoundData)

    nf=m.topology.size[m.topology.dim+1]
    nT=degFs.numB;
    nF=degFv.numB;
    phiT=degFs.phi;
    divphiF=degFv.divphi;
    nPhiT=length(phiT);
    nPhiF=length(divphiF);
    sk=size(kubWeights);
    center=Array{Float64,1}(undef,2);
    nSubCells=compoundData.nSubCells;
    assembledPhiT=compoundData.assembledPhi[degFs.femType];
    assembledPhiF=compoundData.assembledPhi[degFv.femType];
    nCompoundPhiT=compoundData.nCompoundPhi[degFs.femType];
    nCompoundPhiF=compoundData.nCompoundPhi[degFv.femType];
    assembledPhiT=compoundData.assembledPhi[degFs.femType];
    assembledPhiF=compoundData.assembledPhi[degFv.femType];
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end
    mt=m.meshType;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    lS=zeros(length(assembledPhiT), length(assembledPhiF));

    sk=size(kubWeights)
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,nPhiF),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    nquadPhiF=compoundData.nquadPhi[degFv.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights);
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiF_edge=initJacobi((m.geometry.dim,nPhiF),sq);

    for k in 1:nf
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        getSubCells!(subcoord, coord, center, compoundData);
        assemblePhi!(assembledPhiT, subcoord, degFs, m, J, dJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF, subcoord, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        globalNumT=l2g(degFs,k);
        globalNumF=l2g(degFv,k);
        for subCell in 1:nSubCells
            for j in 1:nCompoundPhiF
                for i in 1:nCompoundPhiT
                    currentval=0.0;
                    for subi in 1:nPhiT
                        if assembledPhiT[i][subi,subCell]!=0.0
                            for subj in 1:nPhiF
                                if assembledPhiF[j][subj,subCell]!=0.0
                                    for r in 1:sk[2]
                                        for l in 1:sk[1]
                                            currentval+=assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                                        kubWeights[l,r]*phiT[subi][l,r]*divphiF[subj][l,r];
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if !isequal(currentval,0.0) || (globalNumT[i]==nT && globalNumF[j]==nF)
                        push!(rows,globalNumT[i]);
                        push!(cols,globalNumF[j]);
                        push!(vals,currentval);
                        #ddJ[1] reicht um Vorzeichen der Determinante zu identifizieren
                    end
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

#skalare Größe mit vektorieller Größe
function assembStiffCompound(degFs::degF{1,:H1}, degFv::degF{2,:H1div}, z::Array{Float64,1}, m::mesh, kubWeights::Array{Float64,2}, kubPoints::Array{Float64,2}, compoundData::compoundData)
    nT=degFs.numB;
    nF=degFv.numB;
    phiT=degFs.phi;
    phiF=degFv.phi;
    divphiF=degFv.divphi;
    nPhiT=length(phiT);
    nPhiF=size(phiF,2);
    center=Array{Float64,1}(undef,2);
    nSubCells=compoundData.nSubCells;
    nCompoundPhiT=compoundData.nCompoundPhi[degFs.femType];
    nCompoundPhiF=compoundData.nCompoundPhi[degFv.femType];
    assembledPhiT=compoundData.assembledPhi[degFs.femType];
    assembledPhiF=compoundData.assembledPhi[degFv.femType];
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end
    mt=m.meshType;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    sk=size(kubWeights)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiF=initJacobi((m.geometry.dim,size(phiF,2)),sk);

    nquadPhiF=compoundData.nquadPhi[degFv.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphiF_edge=initJacobi((m.geometry.dim,size(phiF,2)),sq);

    lS=zeros(length(assembledPhiT),length(assembledPhiF));

    for k in 1:m.topology.size[m.topology.dim+1]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        globalNumT=l2g(degFs,k);
        globalNumF=l2g(degFv,k);

        getSubCells!(subcoord, coord, center, compoundData);
        assemblePhi!(assembledPhiT, subcoord, degFs, m, J, ddJ, phiT, kubPoints, kubWeights, compoundData);
        assemblePhi!(assembledPhiF, subcoord, m, divphiF, J_edge, ddJ_edge, jphiF_edge, nquadPhiF, nquadPoints, quadWeights, compoundData);
        for subCell in 1:nSubCells
            jacobi!(J,ddJ,jphiF,kubPoints,phiF,subcoord[subCell],mt);
            for i in 1:nCompoundPhiT
                for j in 1:nCompoundPhiF
                    currentval=0.0;
                    for subi in 1:nPhiT
                        if assembledPhiT[i][subi,subCell]!=0.0
                            for subj in 1:nPhiF
                                if assembledPhiF[j][subj,subCell]!=0.0
                                    for r in 1:sk[2]
                                        for l in 1:sk[1]
                                            zjphi=0.0
                                            for d in 1:m.geometry.dim
                                                zjphi+=z[d]*jphiF[d,subj][l,r]
                                            end
                                            currentval+=assembledPhiT[i][subi,subCell]*assembledPhiF[j][subj,subCell]*
                                                        kubWeights[l,r]*(ddJ[l,r]/abs(ddJ[l,r]))*phiT[subi][l,r]*zjphi;
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if !isequal(currentval,0.0) || (globalNumT[i]==nT && globalNumF[j]==nF)
                        push!(rows,globalNumT[i]);
                        push!(cols,globalNumF[j]);
                        push!(vals,currentval);
                    end
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end




#Gradient einer skalaren Größe mit dem Gradienten einer skalaren Größe (Poisson-Probleme)
function assembStiff!(p::femProblem, comp::Symbol)
    degF=p.degFBoundary;
    m=p.mesh;
    dphiRef=degF[comp].gradphi;

    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    sk=size(kubWeights)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    nDF=degF[comp].numB;
    S=spzeros(nDF,nDF);
    lS=zeros(size(dphiRef,2), size(dphiRef,2));

    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        for j in 1:size(dphiRef,2)
            for i in 1:size(dphiRef,2)
                for k in 1:sk[2]
                    for l in 1:sk[1]
                        lS[i,j]+=kubWeights[l,k]*dJ[l,k]*(dphiRef[1,i][l,k]*dphiRef[1,j][l,k]+dphiRef[2,i][l,k]*dphiRef[2,j][l,k]);
                    end
                end
            end
        end

        globalNum=l2g(degF[comp],k);
        for j in 1:length(globalNum)
            for i in 1:length(globalNum)
                S[globalNum[i],globalNum[j]]+=lS[i,j];
            end
        end
    end
    p.stiffM[comp]=S[1:degF[comp].num,:];
    return nothing;
end
