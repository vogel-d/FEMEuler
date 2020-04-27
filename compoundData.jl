struct compoundData{method}
    nSubCells::Int64;
    nCompoundPhi::Dict{Symbol,Int64};
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
    quadWeights::Array{Float64,1};
    nquadPhi::Dict{Symbol,Array{Array{Array{Float64,1},2},1}};
    nquadPoints::Array{Array{Float64,2},1};
end

function createCompoundData()
    compoundData{:nomethod}(0,Dict(),Dict());
end

function createCompoundData(method::Symbol,femElements::Set{Symbol},m::mesh)
    if method==:HexToKites
        quadPoints, quadWeights=getQuad(9);
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        return compoundData{:HexToKites}(6,nCompoundPhi,initAssembledPhi(6,femElements,m.meshType,nCompoundPhi), quadWeights, nquadPhi, nquadPoints)
    elseif method==:HexToTris
        quadPoints, quadWeights=getQuad(9); #TODO: find better way to handle g (argument of getQuad)
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        return compoundData{:HexToTris}(12,nCompoundPhi,initAssembledPhi(12,femElements,m.meshType,nCompoundPhi), quadWeights, nquadPhi, nquadPoints)
    end
end
