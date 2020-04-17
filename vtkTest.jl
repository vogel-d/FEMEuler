function vtk(m::mesh, degF::degF{1}, sol::Array{Float64,1}, femType::Symbol, filename::String, name::String="Test")

    pts=m.geometry.coordinates;

    Npts=size(pts,2);

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
    fComp=getElementProperties(femType,m.meshType,mx,my);
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    cvtk=zeros(Float64, nf)
    for k in 1:nf
        cLoc=sol[l2g(degF, k)]
        cvtk[k]=dot(fComp,cLoc);
    end
    vtk_cell_data(vtk, cvtk, name)

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function vtk(m::mesh, degF::degF{2}, sol::Array{Float64,1}, femType::Symbol, filename::String, name::String="Test"; printSpherical::Bool=false)

    pts=m.geometry.coordinates;

    Npts=size(pts,2);

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
    fComp=getElementProperties(femType,m.meshType,mx,my);
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    cvtk=zeros(Float64, m.geometry.dim, nf)
    for k in 1:nf
        cLoc=sol[l2g(degF, k)]
        dJ=jacobi!(J,m,k,mx,my,coord);
        fLoc=fComp
        #fLoc=(1/dJ)*J*fComp
        if printSpherical
            xyz=transformation(m,coord,mx,my)
            lon,lat,r=cart2sphere(xyz[1],xyz[2],xyz[3]);
            cvtk[:,k]=velSp(fLoc*cLoc,lon,lat)
        else
            cvtk[:,k]=fLoc*cLoc;
        end
    end
    if size(cvtk,1)==2
        vtk_cell_data(vtk, cvtk[1,:], name*" x")
        vtk_cell_data(vtk, cvtk[2,:], name*" z")
        vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name)
    else
        vtk_cell_data(vtk, cvtk[1,:], name*" x")
        vtk_cell_data(vtk, cvtk[2,:], name*" y")
        vtk_cell_data(vtk, cvtk[3,:], name*" z")
        vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],cvtk[3,:]), name)
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end


function vtk(m::mesh, rx::Int, ry::Int, degF::degF{1}, sol::Array{Float64,1}, femType::Symbol, filename::String, name::String="Test")

    if m.geometry.dim==2
        mf=refineRectMesh(m,rx,ry,:periodic,:constant);
    elseif m.geometry.dim==3
        mf=refineCubedSphere(m,rx)
    end
    pts=mf.geometry.coordinates;
    Npts=size(pts,2);

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
    inc=mf.topology.incidence["20"];
    off=mf.topology.offset["20"];
    nf=mf.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end

    J=Array{Float64,2}(undef,mf.geometry.dim,mf.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,mf.geometry.dim,mf.meshType);
    nfc=m.topology.size[3];

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    fComp=Array{Array{Float64},1}(undef, rx*ry);
    rcoord=Array{Float64,2}(undef,2,rx*ry)
    z=1;
    for j in 1:2:(2*ry-1)
        for i in 1:2:(2*rx-1)
            rcoord[:,z]=[i/(2*rx),j/(2*ry)];
            z+=1;
        end
    end
    for i in 1:size(rcoord,2)
        fComp[i]=getElementProperties(femType,mf.meshType,rcoord[1,i],rcoord[2,i]);
    end
    incm=mf.topology.incidence["CF"]
    offm=mf.topology.offset["CF"]

    cvtk=zeros(Float64, nf)
    cLoc=zeros(Float64, length(degF.phi))
    for k in 1:nfc
        f=incm[offm[k]:offm[k+1]-1];
        cLoc=sol[l2g(degF, k)]
        for i in 1:length(f)
            cvtk[f[i]]=dot(fComp[i],cLoc);
        end
        fill!(cLoc,0.0);
    end
    vtk_cell_data(vtk, cvtk, name)

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function vtk(m::mesh, rx::Int, ry::Int, degF::degF{2}, sol::Array{Float64,1}, femType::Symbol, filename::String, name::String="Test"; printSpherical::Bool=false)

    if m.geometry.dim==2
        mf=refineRectMesh(m,rx,ry,:periodic,:constant);
    elseif m.geometry.dim==3
        mf=refineCubedSphere(m,rx)
    end
    pts=mf.geometry.coordinates;
    Npts=size(pts,2);

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
    inc=mf.topology.incidence["20"];
    off=mf.topology.offset["20"];
    nf=mf.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end

    J=Array{Float64,2}(undef,mf.geometry.dim,mf.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,mf.geometry.dim,mf.meshType);
    nfc=m.topology.size[3];

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    fComp=Array{Array{Float64},1}(undef, rx*ry);
    rcoord=Array{Float64,2}(undef,2,rx*ry)
    z=1;
    for j in 1:2:(2*ry-1)
        for i in 1:2:(2*rx-1)
            rcoord[:,z]=[i/(2*rx),j/(2*ry)];
            z+=1;
        end
    end
    for i in 1:size(rcoord,2)
        fComp[i]=getElementProperties(femType,mf.meshType,rcoord[1,i],rcoord[2,i]);
    end
    incm=mf.topology.incidence["CF"]
    offm=mf.topology.offset["CF"]

    cvtk=zeros(Float64, nf)
    cLoc=zeros(Float64, length(degF.phi))
    for k in 1:nf
        f=incm[offm[k]:offm[k+1]-1];
        cLoc=sol[l2g(degF, k)]
        for i in 1:length(f)
            dJ=jacobi!(J,m,k,rcoord[1,i],rcoord[2,i],coord);
            fLoc=(1/dJ)*J*fComp
            if printSpherical
                xyz=transformation(m,coord,rcoord[1,i],rcoord[2,i])
                lon,lat,r=cart2sphere(xyz[1],xyz[2],xyz[3]);
                cvtk[:,f[i]]=velSp(fLoc*cLoc,lon,lat)
            else
                cvtk[:,f[i]]=fLoc*cLoc;
            end
        end
    end
    if size(cvtk,1)==2
        vtk_cell_data(vtk, cvtk[1,:], name*" x")
        vtk_cell_data(vtk, cvtk[2,:], name*" z")
        vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name)
    else
        vtk_cell_data(vtk, cvtk[1,:], name*" x")
        vtk_cell_data(vtk, cvtk[2,:], name*" y")
        vtk_cell_data(vtk, cvtk[3,:], name*" z")
        vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],cvtk[3,:]), name)
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end


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
