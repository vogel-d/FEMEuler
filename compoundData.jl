struct compoundData
    nSubCells::Int64;
    getSubCells::Array{Array{Array{Float64,2},1},1};
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
end

function createCompoundData()
    compoundData(0,:zero,Dict());
end

function createCompoundData(method::Symbol,femElements::Set{Symbol},m::mesh)
    if method==:HexToKites
        return compoundData(6,getSubCellsHexToKites(m),assemblePhiHexToKites(femElements))
    elseif method==:HexToTris
        return compoundData(12,getSubCellsHexToTris(m),assemblePhiHexToTris(femElements))
    end
end
