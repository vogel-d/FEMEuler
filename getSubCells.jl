function getSubCellsHexToKites!(subcoord,c::Array{Float64,2})
    #beginning from lower right kite
    middle=(1/6).*sum(c,dims=2);
    subcoord[1]=[0.5*(c[1,2]+c[1,3]) c[1,3] 0.5*(c[1,3]+c[1,4]) middle[1];
                 0.5*(c[2,2]+c[2,3]) c[2,3] 0.5*(c[2,3]+c[2,4]) middle[2]];

    subcoord[2]=[0.5*(c[1,3]+c[1,4]) c[1,4] 0.5*(c[1,4]+c[1,5]) middle[1];
                 0.5*(c[2,3]+c[2,4]) c[2,4] 0.5*(c[2,4]+c[2,5]) middle[2]];

    subcoord[3]=[0.5*(c[1,4]+c[1,5]) c[1,5] 0.5*(c[1,5]+c[1,6]) middle[1];
                 0.5*(c[2,4]+c[2,5]) c[2,5] 0.5*(c[2,5]+c[2,6]) middle[2]];

    subcoord[4]=[0.5*(c[1,5]+c[1,6]) c[1,6] 0.5*(c[1,6]+c[1,1]) middle[1];
                 0.5*(c[2,5]+c[2,6]) c[2,6] 0.5*(c[2,6]+c[2,1]) middle[2]];

    subcoord[5]=[0.5*(c[1,6]+c[1,1]) c[1,1] 0.5*(c[1,1]+c[1,2]) middle[1];
                 0.5*(c[2,6]+c[2,1]) c[2,1] 0.5*(c[2,1]+c[2,2]) middle[2]];

    subcoord[6]=[0.5*(c[1,1]+c[1,2]) c[1,2] 0.5*(c[1,2]+c[1,3]) middle[1];
                 0.5*(c[2,1]+c[2,2]) c[2,2] 0.5*(c[2,2]+c[2,3]) middle[2]];
end
