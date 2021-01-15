#Plotten auf verfeinerten Mesh (aktuell nur Rechteck-Mesh)

function refined_vtk(p::femProblem, r::Int, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String; printSpherical::Bool=false)
    m=p.mesh;
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
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/../../VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    sol=p.solution[tend];
    fComp=Array{Array{Float64},1}(undef, r^2);
    rcoord=Array{Float64,2}(undef,2,r^2)
    z=1;
    for j in 1:2:(2*r-1)
        for i in 1:2:(2*r-1)
            rcoord[:,z]=[i/(2*r),j/(2*r)];
            z+=1;
        end
    end

    for l in 1:length(comp)
        solc=getfield(sol,comp[l]);
        for i in 1:size(rcoord,2)
                fComp[i]=getElementProperties(p.femType[comp[l]][1],m.meshType,m.geometry.dim,rcoord[1,i],rcoord[2,i]);
        end
        if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
            cvtk=zeros(Float64, nf*r^2)
            cLoc=zeros(Float64, length(p.degFBoundary[p.femType[comp[l]][1]].phi))
            zk=1
            for k in 1:nf
                f=zk:(zk+r^2-1);
                zk+=r^2
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                for i in 1:length(f)
                    cvtk[f[i]]=dot(fComp[i],cLoc);
                end
                fill!(cLoc,0.0);
            end
            vtk_cell_data(vtk, cvtk, name[l])
        else
            cvtk=zeros(Float64, m.geometry.dim, nf*r^2)
            cLoc=zeros(Float64, size(p.degFBoundary[p.femType[comp[l]][1]].phi,2))
            zk=1
            for k in 1:nf
                f=zk:(zk+r^2-1);
                zk+=r^2
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                for i in 1:length(f)
                    dJ=jacobi!(J,m,k,rcoord[1,i],rcoord[2,i],coord);
                    fLoc=(1/dJ)*J*fComp[i]
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
                vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
                vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
                vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf*r^2)), name[l])
            else
                vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
                vtk_cell_data(vtk, cvtk[2,:], name[l]*" y")
                vtk_cell_data(vtk, cvtk[3,:], name[l]*" z")
                vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],cvtk[3,:]), name[l])
            end
            fill!(cLoc,0.0);
        end
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function refined_vtk(p::femProblem, rx::Int, ry::Int, t::Array{Float64,1}, comp::Array{Symbol,1}, name::Array{String,1}, filename::String; printSpherical::Bool=false)

    m=p.mesh;
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
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/../../VTK/"*filename;

    fComp=Array{Array{Float64},1}(undef, r^2);
    rcoord=Array{Float64,2}(undef,2,r^2)
    z=1;
    for j in 1:2:(2*r-1)
        for i in 1:2:(2*r-1)
            rcoord[:,z]=[i/(2*r),j/(2*r)];
            z+=1;
        end
    end

    outfiles = paraview_collection(vtk_filename_noext) do pvd
        for it in 1:length(t)
            vtk = vtk_grid(@sprintf("%s_%02i", vtk_filename_noext, it), pts, cells)
            sol=p.solution[t[it]];

            for l in 1:length(comp)
                solc=getfield(sol,comp[l]);
                for i in 1:size(rcoord,2)
                        fComp[i]=getElementProperties(p.femType[comp[l]][1],m.meshType,m.geometry.dim,rcoord[1,i],rcoord[2,i]);
                end
                if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
                    cvtk=zeros(Float64, nf*r^2)
                    cLoc=zeros(Float64, length(p.degFBoundary[p.femType[comp[l]][1]].phi))
                    zk=1
                    for k in 1:nf
                        f=zk:(zk+r^2-1);
                        zk+=r^2
                        f=incm[offm[k]:offm[k+1]-1];
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        for i in 1:length(f)
                            cvtk[f[i]]=dot(fComp[i],cLoc);
                        end
                        fill!(cLoc,0.0);
                    end
                    vtk_cell_data(vtk, cvtk, name[l])
                else
                    cvtk=zeros(Float64, m.geometry.dim, nf*r^2)
                    cLoc=zeros(Float64, size(p.degFBoundary[p.femType[comp[l]][1]].phi,2))
                    zk=1
                    for k in 1:nf
                        f=zk:(zk+r^2-1);
                        zk+=r^2
                        f=incm[offm[k]:offm[k+1]-1];
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        for i in 1:length(f)
                            dJ=jacobi!(J,p.mesh,k,rcoord[1,i],rcoord[2,i],coord);
                            fLoc=(1/dJ)*J*fComp[i]
                            if printSpherical
                                xyz=transformation(p.mesh,coord,rcoord[1,i],rcoord[2,i])
                                lon,lat,rad=cart2sphere(xyz[1],xyz[2],xyz[3]);
                                cvtk[:,f[i]]=velSp(fLoc*cLoc,lon,lat)
                            else
                                cvtk[:,f[i]]=fLoc*cLoc;
                            end
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
                    fill!(cLoc,0.0);
                end
            end
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, it)
        end
    end

    return outfiles::Vector{String}
end
