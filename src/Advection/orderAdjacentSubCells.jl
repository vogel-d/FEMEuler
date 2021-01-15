function orderAdjacentSubCells!(adjacentSubCells::Array{Int64,2},subcoord1::Array{Array{Float64,2},1},subcoord2::Array{Array{Float64,2},1},adjacentSubCells1::Array{Int64,1},adjacentSubCells2::Array{Int64,1})
    fill!(adjacentSubCells,0); #just for safety, possibly removeable when working
    z=1;
    for adjacentSubCell1 in adjacentSubCells1
        for adjacentSubCell2 in adjacentSubCells2
            if twocommonpoints(subcoord1[adjacentSubCell1],subcoord2[adjacentSubCell2])
                adjacentSubCells[1,z]=adjacentSubCell1;
                adjacentSubCells[2,z]=adjacentSubCell2;
                z+=1
                break;
            end
        end
    end
    in(0,adjacentSubCells) && @error("ordering not successful.")
end

function twocommonpoints(coord1::Array{Float64,2},coord2::Array{Float64,2},atol::Float64=0.0)
    z=0;
    for i in 1:size(coord1,2)
        for j in 1:size(coord2,2)
            if isapprox(coord1[1,i],coord2[1,j],atol=atol) && isapprox(coord1[2,i],coord2[2,j],atol=atol)
                if z>0
                    return true;
                else
                    z+=1;
                end
            end
        end
    end
    return false
end
