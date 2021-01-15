function assemblePhiPre!(p::femProblem)
    compoundData=p.data.compoundData;
    m=p.mesh;
    femElements=keys(p.degFBoundary);

    femElements=Array{Symbol,1}();
    if p.taskRecovery
        for i in collect(keys(p.femType))
            append!(femElements, p.femType[i])
        end
    else
        for i in collect(keys(p.femType))
            push!(femElements, p.femType[i][1]);
        end
    end

    mt=m.meshType;
    nCells=m.topology.size[m.topology.dim+1];

    center=Array{Float64,1}(undef,2);
    nSubCells=compoundData.nSubCells;
    subcoord=Array{Array{Float64,2},1}(undef,nSubCells);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end

    for femElement in femElements
        nCompoundPhi=compoundData.nCompoundPhi[femElement];
        assembledPhi=compoundData.assembledPhi[femElement];

        if femElement==:DG0
            for Cell in 1:nCells
                coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][Cell]:m.topology.offset["20"][Cell+1]-1]]
                getSubCells!(subcoord, coord, center, compoundData);
                assemblePhi!(assembledPhi, compoundData);
                p.data.compoundData.assembledPhiPre[femElement][Cell]=deepcopy(assembledPhi);
            end
        elseif femElement==:RT0
            quadWeights=compoundData.quadWeights;
            sq=length(quadWeights);
            J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
            ddJ_edge=Array{Float64,1}(undef,sq);
            nquadPhi=compoundData.nquadPhi[:RT0];
            nquadPoints=compoundData.nquadPoints;
            nCompoundPhi=compoundData.nCompoundPhi[:RT0];
            assembledPhi=compoundData.assembledPhi[:RT0];
            phi=p.degFBoundary[:RT0].phi
            divphi=p.degFBoundary[:RT0].divphi
            jphi_edge=initJacobi((m.geometry.dim,size(phi,2)),sq);
            for Cell in 1:nCells
                coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][Cell]:m.topology.offset["20"][Cell+1]-1]]
                getSubCells!(subcoord, coord, center, compoundData);
                assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, quadWeights, compoundData);
                p.data.compoundData.assembledPhiPre[femElement][Cell]=deepcopy(assembledPhi);
            end
        else
            @error("unknown finite element space");
        end
    end
end
