struct compoundData{method}
    nSubCells::Int64;
    nCompoundPhi::Dict{Symbol,Int64};
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
    quadWeights::Array{Float64,1};
    nquadPhi::Dict{Symbol,Array{Array{Array{Float64,1},2},1}};
    nquadPoints::Array{Array{Float64,2},1};
    #boundary[i] is dict describing which subelements have an edge
    #that takes part on edge i of compound element
    boundary::Array{Dict{Int64,Int64},1};
    isEdgePeriodic::Array{Bool,1};
    assembledPhiSafe::Dict{Symbol,Array{Array{Float64,2},1}};
    nVerticesSubElement::Int64;
end

function createCompoundData()
    compoundData{:nomethod}(0,Dict(),Dict(),zeros(1),Dict(),[zeros(1,1)],[Dict(1=>1)],Bool[],Dict(),0);
end

function createCompoundData(method::Symbol,femElements::Set{Symbol},m::mesh)
    if method==:HexToKites
        boundary=Array{Dict{Int64,Int64},1}(undef,6)
        for edge in 1:6
            boundary[edge]=Dict(edge=>2, mod(edge,6)+1=>1)
        end
        quadPoints, quadWeights=getQuad(9);
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        return compoundData{:HexToKites}(6,nCompoundPhi,initAssembledPhi(6,femElements,m.meshType,nCompoundPhi),
                                         quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                         assemblePhiHexToKites(femElements), 4)
    elseif method==:RectToKites
        boundary=Array{Dict{Int64,Int64},1}(undef,4)
        for edge in 1:4
            boundary[edge]=Dict(edge=>2, mod(edge,4)+1=>1)
        end
        quadPoints, quadWeights=getQuad(9);
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>4, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        return compoundData{:RectToKites}(4,nCompoundPhi,initAssembledPhi(4,femElements,m.meshType,nCompoundPhi),
                                         quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                         assemblePhiHexToKites(femElements), 4)
    elseif method==:HexToTris
        boundary=Array{Dict{Int64,Int64},1}(undef,6)
        for edge in 1:6
            boundary[edge]=Dict((edge-1)*2+1=>1, edge*2=>1);
        end
        quadPoints, quadWeights=getQuad(9); #TODO: find better way to handle g (argument of getQuad)
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        return compoundData{:HexToTris}(12,nCompoundPhi,initAssembledPhi(12,femElements,m.meshType,nCompoundPhi),
                                        quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                        assemblePhiHexToKites(femElements), 3)
    end
end
