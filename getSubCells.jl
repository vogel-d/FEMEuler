function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, compoundData::compoundData{:HexToTris})
    middle=(1/6).*sum(coord,dims=2); #allocation

    for i in 1:6
        subcoord[(i-1)*2+1]= [coord[1,1+mod(i+1,6)] 0.5*(coord[1,1+mod(i+1,6)]+coord[1,1+mod(i+2,6)]) middle[1];
                              coord[2,1+mod(i+1,6)] 0.5*(coord[2,1+mod(i+1,6)]+coord[2,1+mod(i+2,6)]) middle[2]];

        subcoord[i*2]=       [0.5*(coord[1,1+mod(i+1,6)]+coord[1,1+mod(i+2,6)]) coord[1,1+mod(i+2,6)] middle[1];
                              0.5*(coord[2,1+mod(i+1,6)]+coord[2,1+mod(i+2,6)]) coord[2,1+mod(i+2,6)] middle[2]];
    end
end


function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, compoundData::compoundData{:HexToKites})
    middle=(1/6).*sum(coord,dims=2);
    for subCell in 1:6
        subcoord[subCell]=[0.5*(coord[1,mod(subCell,6)+1]+coord[1,mod(1+subCell,6)+1]) coord[1,mod(1+subCell,6)+1] 0.5*(coord[1,mod(1+subCell,6)+1]+coord[1,mod(2+subCell,6)+1]) middle[1];
                           0.5*(coord[2,mod(subCell,6)+1]+coord[2,mod(1+subCell,6)+1]) coord[2,mod(1+subCell,6)+1] 0.5*(coord[2,mod(1+subCell,6)+1]+coord[2,mod(2+subCell,6)+1]) middle[2]]
    end
end
