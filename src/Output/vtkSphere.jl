function unstructured_vtkSphere(p::femProblem, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String; printSpherical::Bool=false)

    m=p.mesh;
    Npts=size(m.geometry.coordinates,2);
    pts=Array{Float64,2}(undef,2,Npts)

    for i in 1:Npts
        @views xyz=m.geometry.coordinates[:,i]
        lam,phi,r=cart2sphere(xyz[1],xyz[2],xyz[3]);
        pts[:,i]=[lam, phi];
    end



    if m.meshType==4
        celltype = VTKCellTypes.VTK_QUAD;
        mx=0.5;
        my=0.5;
    elseif m.meshType==3
        celltype = VTKCellTypes.VTK_TRIANGLE;
        mx=0.3333333333333333;
        my=0.3333333333333333;
    else
        error("Ungültiger Mesh-Typ")
    end
    if m.geometry.dim==3
        celltype = VTKCellTypes.VTK_POLYGON;
    end
    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nf=m.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end

    fComp=Array{Array{Float64},1}(undef, length(comp))
    for l in 1:length(comp)
        fComp[l]=getElementProperties(p.femType[comp[l]][1],m.meshType,m.geometry.dim,mx,my);
    end
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = pwd()*"/output/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    sol=p.solution[tend];
    for l in 1:length(comp)
        solc=getfield(sol,comp[l]);
        if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
            cvtk=zeros(Float64, nf)
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                cvtk[k]=dot(fComp[l],cLoc);
            end
            vtk_cell_data(vtk, cvtk, name[l])
        else
            cvtk=zeros(Float64, m.geometry.dim, nf)
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                dJ=jacobi!(J,m,k,mx,my,coord);
                fLoc=(1/dJ)*J*fComp[l]
                if printSpherical
                    xyz=transformation(m,coord,mx,my)
                    lon,lat,r=cart2sphere(xyz[1],xyz[2],xyz[3]);
                    cvtk[:,k]=velSp(fLoc*cLoc,lon,lat)
                else
                    cvtk[:,k]=fLoc*cLoc;
                end
            end
            if size(cvtk,1)==2
                vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
                vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
                vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name[l])
            else
                vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
                vtk_cell_data(vtk, cvtk[2,:], name[l]*" y")
                vtk_cell_data(vtk, cvtk[3,:], name[l]*" z")
                vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],cvtk[3,:]), name[l])
            end
        end
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end
