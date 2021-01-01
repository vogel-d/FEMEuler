function initAssembledPhi(nSubCells::Int64, femElements::Set{Symbol}, meshType::Int64, nCompoundPhi::Dict{Symbol,Int64})
    assembledPhi=Dict();
    dim=2; #compound implementation for 2D only (at the moment)
    for femElement in femElements
        phi,size=getElementProperties(femElement,meshType,dim);
        assembledPhi[femElement]=Array{Array{Float64,2},1}(undef,nCompoundPhi[femElement]);
        for i in 1:nCompoundPhi[femElement]
            assembledPhi[femElement][i]=Array{Float64,2}(undef,size[2],nSubCells)
        end
    end
    return assembledPhi;
end

function initAssembledPhiPre(mesh::mesh, nSubCells::Int64, femElements::Set{Symbol}, nCompoundPhi::Dict{Symbol,Int64})
    meshType=mesh.meshType;
    nCells=mesh.topology.size[mesh.topology.dim+1];
    assembledPhiPre=Dict();
    dim=2; #compound implementation for 2D only (at the moment)

    for femElement in femElements
        phi,size=getElementProperties(femElement,meshType,dim);
        assembledPhiPre[femElement]=Array{Array{Array{Float64,2},1},1}(undef,nCells);
        for Cell in 1:nCells
            assembledPhiPre[femElement][Cell]=Array{Array{Float64,2},1}(undef,nCompoundPhi[femElement])
            for i in 1:nCompoundPhi[femElement]
                assembledPhiPre[femElement][Cell][i]=Array{Float64,2}(undef,size[2],nSubCells)
            end
        end
    end
    return assembledPhiPre;
end
