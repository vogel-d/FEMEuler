function getSubCellsHexToKites(m::mesh)
    subCells=Array{Array{Array{Float64,2},1},1}();
    subCellsCell=[[0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  ]
    for Cell in 1:m.topology.size[m.geometry.dim+1]
        c=m.geometry.coordinates[:,m.topology.incidence["20"][((Cell-1)*6+1):(Cell*6)]];
        middle=(1/6).*sum(c,dims=2);
        for subCell in 1:6
            subCellsCell[subCell]=[0.5*(c[1,mod(subCell,6)+1]+c[1,mod(1+subCell,6)+1]) c[1,mod(1+subCell,6)+1] 0.5*(c[1,mod(1+subCell,6)+1]+c[1,mod(2+subCell,6)+1]) middle[1];
                                   0.5*(c[2,mod(subCell,6)+1]+c[2,mod(1+subCell,6)+1]) c[2,mod(1+subCell,6)+1] 0.5*(c[2,mod(1+subCell,6)+1]+c[2,mod(2+subCell,6)+1]) middle[2]]
        end
        push!(subCells,subCellsCell);
    end
    return subCells;
end


function getSubCellsHexToTris(m::mesh)
    subCells=Array{Array{Array{Float64,2},1},1}();
    subCellsCell=[[0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  [0.0 0.0 0.0 0.0;
                   0.0 0.0 0.0 0.0],
                  ]
    for Cell in 1:m.topology.size[m.geometry.dim+1]
        c=m.geometry.coordinates[:,m.topology.incidence["20"][((Cell-1)*6+1):(Cell*6)]];
        middle=(1/6).*sum(c,dims=2);

        cyclicstartlowerright=[3,4,5,6,1,2];
        for i in 1:6
            current=cyclicstartlowerright[i];
            next=1+mod(current,6);

            subCellsCell[(i-1)*2+1]= [c[1,current] 0.5*(c[1,current]+c[1,next]) middle[1];
                                      c[2,current] 0.5*(c[2,current]+c[2,next]) middle[2]];

            subCellsCell[i*2]=       [0.5*(c[1,current]+c[1,next]) c[1,next] middle[1];
                                      0.5*(c[2,current]+c[2,next]) c[2,next] middle[2]];
        end

        push!(subCells,subCellsCell);
    end
    return subCells;
end
