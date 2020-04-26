struct compoundData{method}
    nSubCells::Int64;
    nCompoundPhi::Dict{Symbol,Int64};
    getSubCells::Array{Array{Array{Float64,2},1},1};
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
end

function createCompoundData()
    compoundData(:nomethod,0,[[zeros(1,1)]],Dict());
end

function createCompoundData(method::Symbol,femElements::Set{Symbol},m::mesh)
    if method==:HexToKites
        return compoundData{:HexToKites}(6,Dict(:RT0=>6, :DG0=>1),getSubCellsHexToKites(m),assemblePhiHexToKites(femElements))
    elseif method==:HexToTris
        return compoundData{:HexToTris}(12,Dict(:RT0=>6, :DG0=>1),getSubCellsHexToTris(m),assemblePhiHexToTris(femElements))
    end
end
