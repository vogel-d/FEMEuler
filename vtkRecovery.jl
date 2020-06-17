function vtkRecovery(m::mesh, r::Int, sol::Array{Float64,1}, femType::Symbol, ::Val{1}, filename::String, name::String="Test")
    if r==2
        refCoordE=[0.5 1.0 0.5 0.0;
                   0.0 0.5 1.0 0.5]

        refCoordF=[0.5; 0.5]

        refInd=[1,5,9,8, 5,2,6,9, 8,9,7,4, 9,6,3,7]
    else
        error("Für das Plotten der $r-fachen Verfeinerung, müssen die Koordinaten und Inzidenzen für das Referenzelement noch angegeben werden.")
    end
    nrCE=size(refCoordE,2)
    nrCF=size(refCoordF,2)
    nf=m.topology.size[3];
    ne=m.topology.size[3];
    nv=m.topology.size[1]
    pts=zeros(Float64,m.geometry.dim,m.topology.size[1]+nf*(nrCE+nrCF))
    pts[:,1:nv]=m.geometry.coordinates;

    m.geometry.dim==3 ? celltype = VTKCellTypes.VTK_POLYGON : celltype = VTKCellTypes.VTK_QUAD;
    mx=0.5; my=0.5;

    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    ze=1; zf=1;
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        coord=@views m.geometry.coordinates[:,inc[off[k]:off[k+1]-1]]
        for i in 1:nrCE
            pts[:,nv+ze]=transformation(m,coord,refCoordE[1,i],refCoordE[2,i])
            push!(inds,nv+ze);
            ze+=1
        end
        for i in 1:nrCF
            pts[:,nv+nf*nrCE+zf]=transformation(m,coord,refCoordF[1,i],refCoordF[2,i])
            push!(inds,nv+nf*nrCE+zf);
            zf+=1
        end
        for i in 1:m.meshType:r^2*m.meshType
            push!(cells, MeshCell(celltype, inds[refInd[i:i+m.meshType-1]]))
        end
    end

    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    fComp=Array{Array{Float64},1}(undef, r^2);
    rcoord=Array{Float64,2}(undef,2,r^2)
    z=1;
    for j in 1:2:(2*r-1)
        for i in 1:2:(2*r-1)
            rcoord[:,z]=[i/(2*r),j/(2*r)];
            z+=1;
        end
    end

    phi=getPhiRecovery(Val(femType))
    for i in 1:r^2
        fComp[i]=similar(phi,Float64);
    end
    nPhi=length(phi)
    cvtk=zeros(Float64, nf*r^2)
    cLoc=zeros(Float64, nPhi)
    zk=1
    for k in 1:nf
        coord=@views m.geometry.coordinates[:,inc[off[k]:off[k+1]-1]]
        mp=transformation(m, coord, mx, my)
        for i in 1:size(rcoord,2)
            mc=transformation(m, coord, rcoord[1,i],rcoord[2,i])
            t1,t2=getTangentialPlane(mp);
            txyz=transformRecoveryCoord(mp,t1,t2,mc)
            for l=1:nPhi
                fComp[i][l]=phi[l](txyz);
            end
            #fComp[i]=getElementProperties(femType,mp,mc);
        end
        f=zk:(zk+r^2-1);
        zk+=r^2
        cLoc=sol[(nPhi*(k-1)+1):(nPhi*k)]
        for i in 1:length(f)
            cvtk[f[i]]=dot(fComp[i],cLoc);
        end
        fill!(cLoc,0.0);
    end
    vtk_cell_data(vtk, cvtk, name)

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function vtkRecovery(m::mesh, r::Int, sol::Array{Float64,1}, femType::Symbol, ::Val{2}, filename::String, name::String="Test"; printSpherical::Bool=false)
    if r==2
        refCoordE=[0.5 1.0 0.5 0.0;
                   0.0 0.5 1.0 0.5]

        refCoordF=[0.5; 0.5]

        refInd=[1,5,9,8, 5,2,6,9, 8,9,7,4, 9,6,3,7]
    else
        error("Für das Plotten der $r-fachen Verfeinerung, müssen die Koordinaten und Inzidenzen für das Referenzelement noch angegeben werden.")
    end
    nrCE=size(refCoordE,2)
    nrCF=size(refCoordF,2)
    nf=m.topology.size[3];
    ne=m.topology.size[3];
    nv=m.topology.size[1]
    pts=zeros(Float64,m.geometry.dim,m.topology.size[1]+nf*(nrCE+nrCF))
    pts[:,1:nv]=m.geometry.coordinates;

    m.geometry.dim==3 ? celltype = VTKCellTypes.VTK_POLYGON : celltype = VTKCellTypes.VTK_QUAD;
    mx=0.5; my=0.5;

    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    ze=1; zf=1;
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        coord=@views m.geometry.coordinates[:,inc[off[k]:off[k+1]-1]]
        for i in 1:nrCE
            pts[:,nv+ze]=transformation(m,coord,refCoordE[1,i],refCoordE[2,i])
            push!(inds,nv+ze);
            ze+=1
        end
        for i in 1:nrCF
            pts[:,nv+nf*nrCE+zf]=transformation(m,coord,refCoordF[1,i],refCoordF[2,i])
            push!(inds,nv+nf*nrCE+zf);
            zf+=1
        end
        for i in 1:m.meshType:r^2*m.meshType
            push!(cells, MeshCell(celltype, inds[refInd[i:i+m.meshType-1]]))
        end
    end

    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    jcoord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    fComp=Array{Array{Float64},1}(undef, r^2);
    rcoord=Array{Float64,2}(undef,2,r^2)
    z=1;
    for j in 1:2:(2*r-1)
        for i in 1:2:(2*r-1)
            rcoord[:,z]=[i/(2*r),j/(2*r)];
            z+=1;
        end
    end

    phi=getPhiRecovery(Val(femType))
    for i in 1:r^2
        fComp[i]=similar(phi,Float64);
    end
    nPhi=size(phi,2)
    lPhi=length(phi)
    cvtk=zeros(Float64, m.geometry.dim,nf*r^2)
    cLoc=zeros(Float64, nPhi)
    zk=1
    for k in 1:nf
        coord=@views m.geometry.coordinates[:,inc[off[k]:off[k+1]-1]]
        mp=transformation(m, coord, mx, my)
        for i in 1:size(rcoord,2)
            mc=transformation(m, coord, rcoord[1,i],rcoord[2,i])
            t1,t2=getTangentialPlane(mp);
            txyz=transformRecoveryCoord(mp,t1,t2,mc)
            for l=1:lPhi
                fComp[i][l]=phi[l](txyz);
            end
            #fComp[i]=getElementProperties(femType,mp,mc);
        end
        f=zk:(zk+r^2-1);
        zk+=r^2
        cLoc=sol[(nPhi*(k-1)+1):(nPhi*k)]
        for i in 1:length(f)
            dJ=jacobi!(J,m,k,rcoord[1,i],rcoord[2,i],jcoord);
            #fLoc=(1/dJ)*J*fComp[i]
            fLoc=fComp[i]
            if printSpherical
                xyz=transformation(m,coord,rcoord[1,i],rcoord[2,i])
                lon,lat,rad=cart2sphere(xyz[1],xyz[2],xyz[3]);
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
