function applyEdge!(phi::Array{Array{Float64,1},2}, ref_phi::Array{Array{Float64,1},2}, eT::Int64, sk::Int64, size::Array{Int64,1}, compoundData::compoundData{:HexToKites})
    if eT>=4
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][1+sk-r]
                end
            end
        end
    else
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][r]
                end
            end
        end
    end
end

function applyEdge!(kubPn::Array{Float64,2}, nquadPoints_ref::Array{Float64,2}, eT::Int64, sk::Int64, compoundData::compoundData{:HexToKites})
    if eT>=4
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[1+sk-r];
        end
    else
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[r];
        end
    end
end

function applyEdge!(phi::Array{Array{Float64,1},2}, ref_phi::Array{Array{Float64,1},2}, eT::Int64, sk::Int64, size::Array{Int64,1}, compoundData::compoundData{:HexToTris})
    if eT>=4
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][1+sk-r]
                end
            end
        end
    else
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][r]
                end
            end
        end
    end
end

function applyEdge!(kubPn::Array{Float64,2}, nquadPoints_ref::Array{Float64,2}, eT::Int64, sk::Int64, compoundData::compoundData{:HexToTris})
    if eT>=4
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[1+sk-r];
        end
    else
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[r];
        end
    end
end

function applyEdge!(phi::Array{Array{Float64,1},2}, ref_phi::Array{Array{Float64,1},2}, eT::Int64, sk::Int64, size::Array{Int64,1}, compoundData::compoundData{:RectToKites})
    if eT>=3
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][1+sk-r]
                end
            end
        end
    else
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][r]
                end
            end
        end
    end
end

function applyEdge!(kubPn::Array{Float64,2}, nquadPoints_ref::Array{Float64,2}, eT::Int64, sk::Int64, compoundData::compoundData{:RectToKites})
    if eT>=3
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[1+sk-r];
        end
    else
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[r];
        end
    end
end

function applyEdge!(phi::Array{Array{Float64,1},2}, ref_phi::Array{Array{Float64,1},2}, eT::Int64, sk::Int64, size::Array{Int64,1}, compoundData::compoundData{:RectToTris})
    if eT>=3
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][1+sk-r]
                end
            end
        end
    else
        for i in 1:size[1]
            for j in 1:size[2]
                for r in 1:sk
                    phi[i,j][r]=ref_phi[i,j][r]
                end
            end
        end
    end
end

function applyEdge!(kubPn::Array{Float64,2}, nquadPoints_ref::Array{Float64,2}, eT::Int64, sk::Int64, compoundData::compoundData{:RectToTris})
    if eT>=3
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[1+sk-r];
        end
    else
        for r in 1:sk
            kubPn[r]=nquadPoints_ref[r];
        end
    end
end
