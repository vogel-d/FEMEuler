function vtkScalar(m::mesh,val::Array{T,1} where T,filename::String)
    coord=m.geometry.coordinates;

    Npts=size(coord,2);
    pts=coord;

    cells = MeshCell[]
    celltype = VTKCellTypes.VTK_POLYGON;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nc=m.topology.size[3];

    for k in 1:nc
        inds=inc[off[k]:off[k+1]-1];
        c = MeshCell(celltype, inds)
        push!(cells, c)
    end
    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    vtk_cell_data(vtk, val, "Werte")

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function vtkVec(m::mesh,val::Array{T,2} where T,filename::String)
    coord=m.geometry.coordinates;

    Npts=size(coord,2);
    pts=coord;

    cells = MeshCell[]
    celltype = VTKCellTypes.VTK_POLYGON;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nc=m.topology.size[3];

    for k in 1:nc
        inds=inc[off[k]:off[k+1]-1];
        c = MeshCell(celltype, inds)
        push!(cells, c)
    end
    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    for i in 1:size(val,1)
        vtk_cell_data(vtk, val[:,i], "Werte $i")
    end
    size(val,1)==2 && vtk_cell_data(vtk, (val[1,:],val[2,:],zeros(Float64,size(val,2))), "Werte")
    size(val,1)==3 && vtk_cell_data(vtk, (val[1,:],val[2,:],val[3,:]), "Werte")

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end
