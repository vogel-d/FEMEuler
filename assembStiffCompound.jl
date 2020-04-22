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
    nSubCells=compoundData.nSubCells;
    assembledPhiT=compoundData.assembledPhi[degFs.femType];
    assembledPhiF=compoundData.assembledPhi[degFv.femType];
    mt=m.meshType;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    # hier ist Assemblieren nur einmal nötig, da sich die Determinanten der Jacobi-Matrix jeweils
    # wegkürzt, weshalb die lokale Steifigkeitsmatrix unabhängig von den Koordinaten der jeweiligen Fläche ist.

    lS=zeros(length(assembledPhiT), length(assembledPhiF));
    for i in 1:length(assembledPhiT)
        for j in 1:length(assembledPhiF)
            currentval=0.0;
            for subCell in 1:nSubCells
                for subi in 1:nPhiT
                    if assembledPhiT[i][subi,subCell]!=0
                        for subj in 1:nPhiF
                            if assembledPhiF[j][subj,subCell]!=0
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
            end
            lS[i,j] = currentval;
        end
    end
    sk=size(kubWeights)
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    checkdet=true;

    for k in 1:nf
        globalNumT=l2g(degFs,k);
        globalNumF=l2g(degFv,k);
        subcoord=compoundData.getSubCells[k];
        for subCell in 1:nSubCells
            jacobi!(J,dJ,kubPoints,subcoord[subCell],mt);
            if (dJ[1]<0 && checkdet)
                @warn("determinant is negative but not considered")
            end
        end
        for j in 1:length(globalNumF)
            for i in 1:length(globalNumT)
                if !isequal(lS[i,j],0.0) || (globalNumT[i]==nT && globalNumF[j]==nF)
                    push!(rows,globalNumT[i]);
                    push!(cols,globalNumF[j]);
                    push!(vals,lS[i,j]);
                    #ddJ[1] reicht um Vorzeichen der Determinante zu identifizieren
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
    nPhiT=length(phiT);
    nPhiF=size(phiF,2);
    nSubCells=compoundData.nSubCells;
    assembledPhiT=compoundData.assembledPhi[degFs.femType];
    assembledPhiF=compoundData.assembledPhi[degFv.femType];
    mt=m.meshType;

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    sk=size(kubWeights)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiF=initJacobi((m.geometry.dim,size(phiF,2)),sk);

    lS=zeros(length(assembledPhiT),length(assembledPhiF));

    for k in 1:m.topology.size[m.topology.dim+1]
        subcoord=compoundData.getSubCells[k];
        for i in 1:length(assembledPhiT)
            for j in 1:length(assembledPhiF)
                currentval=0.0;
                for subCell in 1:nSubCells
                    jacobi!(J,ddJ,jphiF,kubPoints,phiF,subcoord[subCell],mt);
                    for subi in 1:nPhiT
                        if assembledPhiT[i][subi,subCell]!=0
                            for subj in 1:nPhiF
                                if assembledPhiF[j][subj,subCell]!=0
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
                end
                lS[i,j] = currentval;
            end
        end

        globalNumT=l2g(degFs,k);
        globalNumF=l2g(degFv,k);

        for j in 1:length(globalNumF)
            for i in 1:length(globalNumT)
                gi=globalNumT[i];
                gj=globalNumF[j];
                if !isequal(lS[i,j],0.0) || (gi==nT && gj==nF)
                    push!(rows,gi);
                    push!(cols,gj);
                    push!(vals,lS[i,j]);
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
