#Interface zum Assemblieren der globalen Massenmatrix
function assembMassCompound!(p::femProblem)
    degF=p.degFBoundary;
    s=Array{Symbol,1}();
    if p.taskRecovery
        for i in collect(keys(p.femType))
            append!(s, p.femType[i])
        end
    else
        for i in collect(keys(p.femType))
            push!(s, p.femType[i][1]);
        end
    end
    comp=Set{Symbol}(s);
    for i in comp
        @time massM=assembMassCompound(p.degFBoundary[i],p.mesh,p.kubPoints,p.kubWeights,p.compoundData);
        n=p.degFBoundary[i].num
        p.massM[i]=lu(massM[1:n,1:n]);
        p.massMBoundary[i]=lu(massM);
    end
end


#Funktion zum Assemblieren der globalen Massematrix für einen Finite-Elemente-Raum mit skalaren Ansatzfunktionen
function assembMassCompound(degF::degF{1,:H1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, compoundData::compoundData)
    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    phi=degF.phi;
    nPhiSubElement=length(phi);
    sk=size(kubWeights);
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    key="20";
    mt=m.meshType;
    center=Array{Float64,1}(undef,2);
    nSubCells=compoundData.nSubCells;
    nCompoundPhi=compoundData.nCompoundPhi[degF.femType];
    assembledPhi=compoundData.assembledPhi[degF.femType];
    assembledPhiPre=compoundData.assembledPhiPre[degF.femType];
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    for k in 1:m.topology.size[m.topology.dim+1]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        getSubCells!(subcoord, coord, center, compoundData);
        gvertices=l2g(degF,k);
        #assemblePhi!(assembledPhi, compoundData);
        assembledPhi=assembledPhiPre[k];
        for subCell in 1:nSubCells
            jacobi!(J,dJ,kubPoints,subcoord[subCell],mt);
            for j in 1:nCompoundPhi
                for i in 1:nCompoundPhi
                    currentval=0.0;
                    for subj in 1:nPhiSubElement
                        if assembledPhi[j][subj,subCell]!=0.0
                            for subi in 1:nPhiSubElement
                                if assembledPhi[i][subi,subCell]!=0.0
                                    for r in 1:sk[2]
                                        for l in 1:sk[1]
                                            currentval+=assembledPhi[i][subi,subCell]*assembledPhi[j][subj,subCell]*
                                                        kubWeights[l,r]*phi[subi][l,r]*phi[subj][l,r]*abs(dJ[l,r]);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if !isequal(currentval,0.0)
                        push!(rows,gvertices[i]);
                        push!(cols,gvertices[j]);
                        push!(vals,currentval);
                    end
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

#Funktion zum Assemblieren der globalen Massematrix für einen Finite-Elemente-Raum mit vektoriellen Ansatzfunktionen
function assembMassCompound(degF::degF{2,:H1div}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, compoundData::compoundData)
#function assembMass(degF::degF{2,S} where S, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    phi=degF.phi;
    divphi=degF.divphi;
    nPhiSubElement=size(phi,2);
    sk=size(kubWeights)
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,nPhiSubElement),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,diff(m.topology.offset["20"][1:2])[1]);
    key="20";
    mt=m.meshType;
    center=Array{Float64,1}(undef,2);
    nSubCells=compoundData.nSubCells;
    nquadPhi=compoundData.nquadPhi[degF.femType];
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    nCompoundPhi=compoundData.nCompoundPhi[degF.femType];
    assembledPhi=compoundData.assembledPhi[degF.femType];
    assembledPhiPre=compoundData.assembledPhiPre[degF.femType];
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    sq=length(quadWeights)
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphi_edge=initJacobi((m.geometry.dim,size(phi,2)),sq);

    for k in 1:m.topology.size[m.topology.dim+1]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]

        getSubCells!(subcoord, coord, center, compoundData);
        gvertices=l2g(degF,k);
        #assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, quadWeights, compoundData);
        assembledPhi=assembledPhiPre[k];
        for subCell in 1:nSubCells
            jacobi!(J,ddJ,jphi,kubPoints,phi,subcoord[subCell],mt);
            for j in 1:nCompoundPhi
                for i in 1:nCompoundPhi
                    currentval=0.0;
                    for subj in 1:nPhiSubElement
                        if assembledPhi[j][subj,subCell]!=0.0
                            for subi in 1:nPhiSubElement
                                if assembledPhi[i][subi,subCell]!=0.0
                                    for r in 1:sk[2]
                                        for l in 1:sk[1]
                                            dotp=0.0;
                                            for d in 1:m.geometry.dim
                                                dotp+=jphi[d,subi][l,r]*jphi[d,subj][l,r];
                                            end
                                            currentval+=assembledPhi[i][subi,subCell]*assembledPhi[j][subj,subCell]*
                                                        kubWeights[l,r]*abs(ddJ[l,r])*dotp
                                        end
                                    end
                                end
                            end
                        end
                    end
                    if !isequal(currentval,0.0)
                        push!(rows,gvertices[i]);
                        push!(cols,gvertices[j]);
                        push!(vals,currentval);
                    end
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end


function assembMassCompound(degF::degF{2,:H1xH1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, compoundData::compoundData)
    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    phi=degF.phi;
    iter=size(phi,2);
    sk=size(kubWeights)
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim,iter),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);
        gvertices=l2g(degF,k);
        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        vecdot=0.0;
                        for d in 1:m.geometry.dim
                            vecdot+=phi[d,i][l,r]*phi[d,j][l,r];
                        end
                        currentval+=kubWeights[l,r]*abs(dJ[l,r])*vecdot
                    end
                end

                if !isequal(currentval,0.0)
                    push!(rows,gvertices[i]);
                    push!(cols,gvertices[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end
