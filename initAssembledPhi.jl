function initAssembledPhi(nSubCells::Int64, femElements::Set{Symbol}, meshType::Int64, nCompoundPhi::Dict{Symbol,Int64})
    assembledPhi=Dict();
    for femElement in femElements
        phi,size=getElementProperties(femElement,meshType);
        assembledPhi[femElement]=Array{Array{Float64,2},1}(undef,nCompoundPhi[femElement]);
        for i in 1:nCompoundPhi[femElement]
            assembledPhi[femElement][i]=Array{Float64,2}(undef,size[2],nSubCells)
        end
    end
    return assembledPhi;
end
