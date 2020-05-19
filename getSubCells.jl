function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, center::Array{Float64,1}, compoundData::compoundData{:HexToTris})
    #middle=(1/6).*sum(coord,dims=2); #allocation
    fill!(center,0.0);
    for i in 1:6
        for j in 1:2
            center[j]+=(1/6)*coord[j,i];
        end
    end

    for i in 1:6
        #subcoord[(i-1)*2+1]= [coord[1,i] 0.5*(coord[1,i]+coord[1,1+mod(i,6)]) middle[1];
        #                      coord[2,i] 0.5*(coord[2,i]+coord[2,1+mod(i,6)]) middle[2]];

        #subcoord[i*2]=       [0.5*(coord[1,i]+coord[1,1+mod(i,6)]) coord[1,1+mod(i,6)] middle[1];
        #                      0.5*(coord[2,i]+coord[2,1+mod(i,6)]) coord[2,1+mod(i,6)] middle[2]];

        for j in 1:2
          subcoord[(i-1)*2+1][j,1]=coord[j,i]
        end

        for j in 1:2
          subcoord[(i-1)*2+1][j,2]=0.5*(coord[j,i]+coord[j,1+mod(i,6)])
        end

        for j in 1:2
          subcoord[(i-1)*2+1][j,3]=center[j]
        end

        for j in 1:2
          subcoord[i*2][j,1]=0.5*(coord[j,i]+coord[j,1+mod(i,6)])
        end

        for j in 1:2
          subcoord[i*2][j,2]=coord[j,1+mod(i,6)]
        end

        for j in 1:2
          subcoord[i*2][j,3]=center[j]
        end
    end
end


function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, center::Array{Float64,1}, compoundData::compoundData{:HexToKites})
    #middle=(1/6).*sum(coord,dims=2);
    fill!(center,0.0);
    for i in 1:6
        for j in 1:2
            center[j]+=(1/6)*coord[j,i];
        end
    end
    for i in 1:6
        for j in 1:2
            subcoord[i][j,1]=0.5*(coord[j,mod(i+4,6)+1]+coord[j,i])
        end

        for j in 1:2
            subcoord[i][j,2]=coord[j,i];
        end

        for j in 1:2
            subcoord[i][j,3]=0.5*(coord[j,i]+coord[j,mod(i,6)+1])
        end

        for j in 1:2
            subcoord[i][j,4]=center[j];
        end
        #subcoord[i]=[0.5*(coord[1,mod(i+4,6)+1]+coord[1,i]) coord[1,i] 0.5*(coord[1,i]+coord[1,mod(i,6)+1]) center[1];
        #             0.5*(coord[2,mod(i+4,6)+1]+coord[2,i]) coord[2,i] 0.5*(coord[2,i]+coord[2,mod(i,6)+1]) center[2]]
    end
end

function getSubCells!(subcoord::Array{Array{Float64,2},1}, coord::Array{Float64,2}, center::Array{Float64,1}, compoundData::compoundData{:RectToKites})
    #middle=(1/6).*sum(coord,dims=2);
    fill!(center,0.0);
    for i in 1:4
        for j in 1:2
            center[j]+=(1/4)*coord[j,i];
        end
    end
    for i in 1:4
        for j in 1:2
            subcoord[i][j,1]=0.5*(coord[j,mod(i+2,4)+1]+coord[j,i])
        end

        for j in 1:2
            subcoord[i][j,2]=coord[j,i];
        end

        for j in 1:2
            subcoord[i][j,3]=0.5*(coord[j,i]+coord[j,mod(i,4)+1])
        end

        for j in 1:2
            subcoord[i][j,4]=center[j];
        end
        #subcoord[i]=[0.5*(coord[1,mod(i+4,6)+1]+coord[1,i]) coord[1,i] 0.5*(coord[1,i]+coord[1,mod(i,6)+1]) center[1];
        #             0.5*(coord[2,mod(i+4,6)+1]+coord[2,i]) coord[2,i] 0.5*(coord[2,i]+coord[2,mod(i,6)+1]) center[2]]
    end
end
