struct compoundData{method}
    nSubCells::Int64;
    nCompoundPhi::Dict{Symbol,Int64};
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
end

function createCompoundData()
    compoundData{:nomethod}(0,Dict(),Dict());
end

function createCompoundData(method::Symbol,femElements::Set{Symbol},m::mesh)
    if method==:HexToKites
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        return compoundData{:HexToKites}(6,nCompoundPhi,initAssembledPhi(6,femElements,m.meshType,nCompoundPhi))
    elseif method==:HexToTris
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        return compoundData{:HexToTris}(12,nCompoundPhi,initAssembledPhi(12,femElements,m.meshType,nCompoundPhi))
    end
end
