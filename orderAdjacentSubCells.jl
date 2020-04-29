function orderAdjacentSubCells!(adjacentSubCells::Array{Int64,2},m::mesh,subcoord1,subcoord2,adjacentSubCells1,adjacentSubCells2)
    coordCell1= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc1]:m.topology.offset["20"][inc1+1]-1]]
    coordCell2= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][inc2]:m.topology.offset["20"][inc2+1]-1]]
    getSubCells!(subcoord1, coordCell1, compoundData);
    getSubCells!(subcoord2, coordCell2, compoundData);

    z=0;
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
end

function twocommonpoints(coord1::Array{Float64,2},coord2::Array{Float64,2})
    commonpoints=Int64[];
    for i in 1:size(coord1,2)
        for commonpoint in findall(coord1[:,i],coord2)
            push!(commonpoints,commonpoint);
        end
    end
    return length(commonpoints)==2
end
