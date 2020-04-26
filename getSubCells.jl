function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, compoundData::compoundData{:HexToTris})
    middle=(1/6).*sum(coord,dims=2); #allocation

    for i in 1:6
        subcoord[(i-1)*2+1]= [coord[1,i] 0.5*(coord[1,i]+coord[1,1+mod(i,6)]) middle[1];
                              coord[2,i] 0.5*(coord[2,i]+coord[2,1+mod(i,6)]) middle[2]];

        subcoord[i*2]=       [0.5*(coord[1,i]+coord[1,1+mod(i,6)]) coord[1,1+mod(i,6)] middle[1];
                              0.5*(coord[2,i]+coord[2,1+mod(i,6)]) coord[2,1+mod(i,6)] middle[2]];
    end
end


function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, compoundData::compoundData{:HexToKites})
    middle=(1/6).*sum(coord,dims=2);
    for i in 1:6
        subcoord[i]=[0.5*(coord[1,mod(i+4,6)+1]+coord[1,i]) coord[1,i] 0.5*(coord[1,i]+coord[1,mod(i,6)+1]) middle[1];
                     0.5*(coord[2,mod(i+4,6)+1]+coord[2,i]) coord[2,i] 0.5*(coord[2,i]+coord[2,mod(i,6)+1]) middle[2]]
    end
end
