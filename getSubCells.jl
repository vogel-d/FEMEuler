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
        for subCell in 1:2:11
            subCellsCell[subCell]=[0.5*(c[1,mod(subCell,6)+1]+c[1,mod(1+subCell,6)+1]) c[1,mod(1+subCell,6)+1] middle[1];
                                   0.5*(c[2,mod(subCell,6)+1]+c[2,mod(1+subCell,6)+1]) c[2,mod(1+subCell,6)+1] middle[2]]

            subCellsCell[subCell+1]=[c[1,mod(1+subCell,6)+1] 0.5*(c[1,mod(1+subCell,6)+1]+c[1,mod(2+subCell,6)+1]) middle[1];
                                     c[2,mod(1+subCell,6)+1] 0.5*(c[2,mod(1+subCell,6)+1]+c[2,mod(2+subCell,6)+1]) middle[2]]

        end
        push!(subCells,subCellsCell);
    end
    return subCells;
end
