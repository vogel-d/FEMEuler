using WriteVTK
using Printf

function unstructured_vtk(p::femProblem, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    m=p.mesh;
    pts=m.geometry.coordinates;

    Npts=size(pts,2);

    m.meshType==4 ? celltype = VTKCellTypes.VTK_QUAD : celltype = VTKCellTypes.VTK_TRIANGLE;
    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nf=m.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end
    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    sol=p.solution[tend];
    for l in 1:length(comp)
        solc=getfield(sol,comp[l]);
        c=p.degFBoundary[p.femType[comp[l]][1]].components;
        if c[1]==0
            cvtk=Float64[mean(solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]) for k in 1:nf]
            vtk_cell_data(vtk, cvtk, name[l])
        else
            xi=findall(c.==1)
            yi=findall(c.==2)
            cvtk1=Float64[mean(solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)[xi]]) for k in 1:nf]
            cvtk2=Float64[mean(solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)[yi]]) for k in 1:nf]
            cvtk3=zeros(Float64,nf)
            vtk_cell_data(vtk, cvtk1, name[l]*" x")
            vtk_cell_data(vtk, cvtk2, name[l]*" z")
            vtk_cell_data(vtk, (cvtk1,cvtk2,cvtk3), name[l])
        end
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function unstructured_vtk(p::femProblem, t::Array{Float64,1}, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    m=p.mesh;
    pts=m.geometry.coordinates;

    Npts=size(pts,2);

    m.meshType==4 ? celltype = VTKCellTypes.VTK_QUAD : celltype = VTKCellTypes.VTK_TRIANGLE;
    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nf=m.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;

    outfiles = paraview_collection(vtk_filename_noext) do pvd
        for it in 1:length(t)
            vtk = vtk_grid(@sprintf("%s_%02i", vtk_filename_noext, it), pts, cells)
            sol=p.solution[t[it]];
            for l in 1:length(comp)
                solc=getfield(sol,comp[l]);
                c=p.degFBoundary[p.femType[comp[l]][1]].components;
                if c[1]==0
                    cvtk=Float64[mean(solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]) for k in 1:nf]
                    vtk_cell_data(vtk, cvtk, name[l])
                else
                    xi=findall(c.==1)
                    yi=findall(c.==2)
                    cvtk1=Float64[mean(solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)[xi]]) for k in 1:nf]
                    cvtk2=Float64[mean(solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)[yi]]) for k in 1:nf]
                    cvtk3=zeros(Float64,nf)
                    vtk_cell_data(vtk, cvtk1, name[l]*" x")
                    vtk_cell_data(vtk, cvtk2, name[l]*" z")
                    vtk_cell_data(vtk, (cvtk1,cvtk2,cvtk3), name[l])
                end
            end
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, it)
        end
    end

    return outfiles::Vector{String}
end
