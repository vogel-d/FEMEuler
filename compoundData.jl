struct compoundData{method}
    nSubCells::Int64;
    nCompoundPhi::Dict{Symbol,Int64};
    assembledPhi::Dict{Symbol,Array{Array{Float64,2},1}};
    quadWeights::Array{Float64,1};
    nquadPhi::Dict{Symbol,Array{Array{Array{Float64,1},2},1}};
    nquadPoints::Array{Array{Float64,2},1};
    #boundary[i] is dict describing which subelements have an edge
    #that takes part on edge i of the compound element
    boundary::Array{Dict{Int64,Int64},1};
    isEdgePeriodic::Array{Bool,1};
    assembledPhiPre::Dict{Symbol,Array{Array{Array{Float64,2},1},1}};
    nVerticesSubElement::Int64;
end

function createCompoundData()
    compoundData{:nomethod}(0,Dict(),Dict(),zeros(1),Dict(),[zeros(1,1)],[Dict(1=>1)],Bool[],Dict(),0);
end

function createCompoundData(method::Symbol,femElements::Set{Symbol},m::mesh)
    if method==:HexToKites
        m.meshType!=4 && error("meshtype is not suitable for quadrilateral subelements")
        boundary=Array{Dict{Int64,Int64},1}(undef,6)
        for edge in 1:6
            boundary[edge]=Dict(edge=>2, mod(edge,6)+1=>1)
        end
        quadPoints, quadWeights=getQuad(9);
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        nVerticesSubElement=4
        nSubCells=6
        return compoundData{:HexToKites}(nSubCells,nCompoundPhi,initAssembledPhi(nSubCells,femElements,m.meshType,nCompoundPhi),
                                         quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                         initAssembledPhiPre(m,nSubCells,femElements,nCompoundPhi), nVerticesSubElement)
    elseif method==:RectToKites
        m.meshType!=4 && error("meshtype is not suitable for quadrilateral subelements")
        boundary=Array{Dict{Int64,Int64},1}(undef,4)
        for edge in 1:4
            boundary[edge]=Dict(edge=>2, mod(edge,4)+1=>1)
        end
        quadPoints, quadWeights=getQuad(9);
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>4, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        nVerticesSubElement=4
        nSubCells=4
        return compoundData{:RectToKites}(nSubCells,nCompoundPhi,initAssembledPhi(nSubCells,femElements,m.meshType,nCompoundPhi),
                                         quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                         initAssembledPhiPre(m,nSubCells,femElements,nCompoundPhi), nVerticesSubElement)
    elseif method==:HexToTris
        m.meshType!=3 && error("meshtype is not suitable for triangular subelements")
        boundary=Array{Dict{Int64,Int64},1}(undef,6)
        for edge in 1:6
            boundary[edge]=Dict((edge-1)*2+1=>1, edge*2=>1);
        end
        quadPoints, quadWeights=getQuad(9); #TODO: find better way to handle g (argument of getQuad)
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>6, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        nVerticesSubElement=3
        nSubCells=12
        return compoundData{:HexToTris}(nSubCells,nCompoundPhi,initAssembledPhi(nSubCells,femElements,m.meshType,nCompoundPhi),
                                        quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                        initAssembledPhiPre(m,nSubCells,femElements,nCompoundPhi), nVerticesSubElement)
    elseif method==:RectToTris
        m.meshType!=3 && error("meshtype is not suitable for triangular subelements")
        boundary=Array{Dict{Int64,Int64},1}(undef,4)
        for edge in 1:4
            boundary[edge]=Dict((edge-1)*2+1=>1, edge*2=>1);
        end
        quadPoints, quadWeights=getQuad(9);
        nquadPhi, nquadPoints=coordTrans(m.meshType, m.normals, [:RT0], 5);
        nCompoundPhi=Dict(:RT0=>4, :DG0=>1);
        femElements=Set([:RT0,:DG0]);
        nVerticesSubElement=3
        nSubCells=8
        return compoundData{:RectToTris}(nSubCells,nCompoundPhi,initAssembledPhi(nSubCells,femElements,m.meshType,nCompoundPhi),
                                         quadWeights, nquadPhi, nquadPoints, boundary, Bool[],
                                         initAssembledPhiPre(m,nSubCells,femElements,nCompoundPhi), nVerticesSubElement)
    else
        error("method of compound finite elements not available")
    end
end
