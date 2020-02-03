using WriteVTK
using Printf

function unstructured_vtk(p::femProblem, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    m=p.mesh;
    pts=m.geometry.coordinates;

    Npts=size(pts,2);

    #m.meshType==4 ? celltype = VTKCellTypes.VTK_QUAD : celltype = VTKCellTypes.VTK_TRIANGLE;
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
        fComp[l]=getElementProperties(p.femType[comp[l]][1],m.meshType,mx,my);
    end
    J=Array{Float64,2}(undef,2,2);
    dJ=0.0;
    coord=Array{Float64,2}(undef,2,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
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
            cvtk=zeros(Float64, 2, nf)
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                dJ=jacobi!(J,m,k,mx,my,coord);
                fLoc=(1/dJ)*J*fComp[l]
                cvtk[:,k]=fLoc*cLoc;
            end
            vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
            vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
            vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name[l])
        end
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function unstructured_vtk(p::femProblem, t::Array{Float64,1}, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    m=p.mesh;
    pts=m.geometry.coordinates;

    Npts=size(pts,2);

    #m.meshType==4 ? celltype = VTKCellTypes.VTK_QUAD : celltype = VTKCellTypes.VTK_TRIANGLE;
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
        fComp[l]=getElementProperties(p.femType[comp[l]][1],m.meshType,mx,my);
    end
    J=Array{Float64,2}(undef,2,2);
    dJ=0.0;
    coord=Array{Float64,2}(undef,2,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;

    outfiles = paraview_collection(vtk_filename_noext) do pvd
        for it in 1:length(t)
            vtk = vtk_grid(@sprintf("%s_%02i", vtk_filename_noext, it), pts, cells)
            sol=p.solution[t[it]];
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
                    cvtk=zeros(Float64, 2, nf)
                    for k in 1:nf
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        dJ=jacobi!(J,m,k,mx,my,coord);
                        fLoc=(1/dJ)*J*fComp[l]
                        cvtk[:,k]=fLoc*cLoc;
                    end
                    vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
                    vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
                    vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name[l])
                end
            end
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, it)
        end
    end

    return outfiles::Vector{String}
end


function unstructured_vtk(p::femProblem, rx::Int, ry::Int, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    mf=refineRectMesh(p.mesh,rx,ry,:periodic,:constant);
    pts=mf.geometry.coordinates;
    Npts=size(pts,2);

    mf.meshType==4 ? celltype = VTKCellTypes.VTK_QUAD : celltype = VTKCellTypes.VTK_TRIANGLE;
    cells = MeshCell[]
    inc=mf.topology.incidence["20"];
    off=mf.topology.offset["20"];
    nf=mf.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end

    J=Array{Float64,2}(undef,2,2);
    dJ=0.0;
    coord=Array{Float64,2}(undef,2,mf.meshType);
    nfc=p.mesh.topology.size[3];

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    sol=p.solution[tend];
    fComp=Array{Array{Float64},1}(undef, rx*ry);
    rcoord=Array{Float64,2}(undef,2,rx*ry)
    z=1;
    for j in 1:2:(2*ry-1)
        for i in 1:2:(2*rx-1)
            rcoord[:,z]=[i/(2*rx),j/(2*ry)];
            z+=1;
        end
    end
    incm=mf.topology.incidence["CF"]
    offm=mf.topology.offset["CF"]
    for l in 1:length(comp)
        solc=getfield(sol,comp[l]);
        for i in 1:size(rcoord,2)
                fComp[i]=getElementProperties(p.femType[comp[l]][1],mf.meshType,rcoord[1,i],rcoord[2,i]);
        end
        if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
            cvtk=zeros(Float64, nf)
            cLoc=zeros(Float64, length(p.degFBoundary[p.femType[comp[l]][1]].phi))
            for k in 1:nfc
                f=incm[offm[k]:offm[k+1]-1];
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                for i in 1:length(f)
                    cvtk[f[i]]=dot(fComp[i],cLoc);
                end
                fill!(cLoc,0.0);
            end
            vtk_cell_data(vtk, cvtk, name[l])
        else
            cvtk=zeros(Float64, 2, nf)
            cLoc=zeros(Float64, size(p.degFBoundary[p.femType[comp[l]][1]].phi,2))
            for k in 1:nfc
                f=incm[offm[k]:offm[k+1]-1];
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                for i in 1:length(f)
                    dJ=jacobi!(J,p.mesh,k,rcoord[1,i],rcoord[2,i],coord);
                    fLoc=(1/dJ)*J*fComp[i]
                    cvtk[:,f[i]]=fLoc*cLoc;
                end
            end
            vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
            vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
            vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name[l])
            fill!(cLoc,0.0);
        end
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function unstructured_vtk(p::femProblem, rx::Int, ry::Int, t::Array{Float64,1}, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    mf=refineRectMesh(p.mesh,rx,ry,:periodic,:constant);
    pts=mf.geometry.coordinates;
    Npts=size(pts,2);

    mf.meshType==4 ? celltype = VTKCellTypes.VTK_QUAD : celltype = VTKCellTypes.VTK_TRIANGLE;
    cells = MeshCell[]
    inc=mf.topology.incidence["20"];
    off=mf.topology.offset["20"];
    nf=mf.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end

    J=Array{Float64,2}(undef,2,2);
    dJ=0.0;
    coord=Array{Float64,2}(undef,2,mf.meshType);
    nfc=p.mesh.topology.size[3];

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;

    fComp=Array{Array{Float64},1}(undef, rx*ry);
    rcoord=Array{Float64,2}(undef,2,rx*ry)

    outfiles = paraview_collection(vtk_filename_noext) do pvd
        for it in 1:length(t)
            vtk = vtk_grid(@sprintf("%s_%02i", vtk_filename_noext, it), pts, cells)
            sol=p.solution[t[it]];

            for l in 1:length(comp)
                solc=getfield(sol,comp[l]);
                for i in 1:size(rcoord,2)
                        fComp[i]=getElementProperties(p.femType[comp[l]][1],mf.meshType,rcoord[1,i],rcoord[2,i]);
                end
                if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
                    cvtk=zeros(Float64, nf)
                    cLoc=zeros(Float64, length(p.degFBoundary[p.femType[comp[l]][1]].phi))
                    for k in 1:nfc
                        f=incm[offm[k]:offm[k+1]-1];
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        for i in 1:length(f)
                            cvtk[f[i]]=dot(fComp[i],cLoc);
                        end
                        fill!(cLoc,0.0);
                    end
                    vtk_cell_data(vtk, cvtk, name[l])
                else
                    cvtk=zeros(Float64, 2, nf)
                    cLoc=zeros(Float64, size(p.degFBoundary[p.femType[comp[l]][1]].phi,2))
                    for k in 1:nfc
                        f=incm[offm[k]:offm[k+1]-1];
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        for i in 1:length(f)
                            dJ=jacobi!(J,p.mesh,k,rcoord[1,i],rcoord[2,i],coord);
                            fLoc=(1/dJ)*J*fComp[i]
                            cvtk[:,f[i]]=fLoc*cLoc;
                        end
                    end
                    vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
                    vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
                    vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name[l])
                    fill!(cLoc,0.0);
                end
            end
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, it)
        end
    end

    return outfiles::Vector{String}
end
