using WriteVTK
using Printf

function compound_unstructured_vtk(p::femProblem, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    m=p.mesh;
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
    if typeof(p.compoundData).parameters[1]!=:nomethod
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
        fComp[l]=getCompoundElementProperties(p.femType[comp[l]][1],m.meshType,p.compoundData.assembledPhi[p.femType[comp[l]][1]],mx,my,p.compoundData);
    end
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    sol=p.solution[tend];
    for l in 1:length(comp)
        solc=getfield(sol,comp[l]);
        if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
            cvtk=zeros(Float64, nf)
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                #cvtk[k]=dot(fComp[l],cLoc);
                cvtk[k]=mean(cLoc);
            end
            vtk_cell_data(vtk, cvtk, name[l])
        else
            cvtk=zeros(Float64, m.geometry.dim, nf)
            fCompoundLoc=zeros(Float64, m.geometry.dim, p.compoundData.nCompoundPhi[p.femType[comp[l]][1]])
            center=Array{Float64,1}(undef,2);
            subcoord=Array{Array{Float64,2},1}(undef,p.compoundData.nSubCells);
            assembledPhi=p.compoundData.assembledPhi[p.femType[comp[l]][1]]
            divphi=p.degFBoundary[p.femType[comp[l]][1]].divphi;
            sq=length(p.compoundData.quadWeights)
            J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
            ddJ_edge=Array{Float64,1}(undef,sq);
            jphi_edge=initJacobi((m.geometry.dim,size(phi,2)),sq);
            for i in 1:nSubCells
                #fill! causes mutating all entries of subcoord when changing a single entry
                subcoord[i]=Array{Float64,2}(undef,2,p.compoundData.nVerticesSubElement);
            end
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
                getSubCells!(subcoord, coord, center, p.compoundData);
                assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, p.quadWeights, p.compoundData);
                fill!(fCompoundLoc,0.0)
                for subCell in 1:p.compoundData.nSubCells
                    dJ=jacobi!(J,mx,my,subcoord[subCell],mt);
                    fLoc=(1/dJ)*J*fComp[l]
                    for i in 1:length(assembledPhi)
                        for subi in 1:size(phi,2)
                            fCompoundLoc[:,i]+=(1/p.compoundData.nSubCells)*(assembledPhi[i][subi,subCell]*fLoc[:,subi]);
                        end
                    end
                end
                cvtk[:,k]=fCompoundLoc*cLoc;
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

function compound_unstructured_vtk(p::femProblem, t::Array{Float64,1}, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    m=p.mesh;
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
    if typeof(p.compoundData).parameters[1]!=:nomethod
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
        fComp[l]=getCompoundElementProperties(p.femType[comp[l]][1],m.meshType,mx,my,p.compoundData);
    end
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

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
                        #cvtk[k]=dot(fComp[l],cLoc);
                        cvtk[k]=mean(cLoc);
                    end
                    vtk_cell_data(vtk, cvtk, name[l])
                else
                    #=
                    cvtk=zeros(Float64, m.geometry.dim, nf)
                    for k in 1:nf
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        #dJ=jacobi!(J,m,k,mx,my,coord);
                        #fLoc=(1/dJ)*J*fComp[l]
                        #cvtk[:,k]=fLoc*cLoc;
                        cvtk[:,k]=[mean(cLoc), mean(cLoc)];
                    end
                    =#
                    cvtk=zeros(Float64, m.geometry.dim, nf)
                    fCompoundLoc=zeros(Float64, m.geometry.dim, p.compoundData.nCompoundPhi[p.femType[comp[l]][1]])
                    center=Array{Float64,1}(undef,2);
                    subcoord=Array{Array{Float64,2},1}(undef,p.compoundData.nSubCells);
                    assembledPhi=p.compoundData.assembledPhi[p.femType[comp[l]][1]]
                    divphi=p.degFBoundary[p.femType[comp[l]][1]].divphi;
                    sq=length(p.compoundData.quadWeights);
                    nquadPhi=p.compoundData.nquadPhi[p.femType[comp[l]][1]];
                    nSubCells=p.compoundData.nSubCells;
                    nquadPoints=p.compoundData.nquadPoints;
                    quadWeights=p.compoundData.quadWeights;
                    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
                    ddJ_edge=Array{Float64,1}(undef,sq);
                    jphi_edge=initJacobi((m.geometry.dim,size(fComp[l],2)),sq);
                    for i in 1:nSubCells
                        #fill! causes mutating all entries of subcoord when changing a single entry
                        subcoord[i]=Array{Float64,2}(undef,2,p.compoundData.nVerticesSubElement);
                    end
                    for k in 1:nf
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
                        getSubCells!(subcoord, coord, center, p.compoundData);
                        assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, quadWeights, p.compoundData);
                        fill!(fCompoundLoc,0.0)
                        for subCell in 1:nSubCells
                            dJ=jacobi!(J,mx,my,subcoord[subCell],m.meshType);
                            fLoc=(1/dJ)*J*fComp[l]
                            for i in 1:length(assembledPhi)
                                for subi in 1:size(fComp[l],2)
                                    fCompoundLoc[:,i]+=(1/nSubCells)*(assembledPhi[i][subi,subCell]*fLoc[:,subi]);
                                end
                            end
                        end
                        cvtk[:,k]=fCompoundLoc*cLoc;
                        #println(fCompoundLoc)
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
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, it)
        end
    end

    return outfiles::Vector{String}
end

function split_unstructured_vtk(p::femProblem, tend::Float64, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    compound_m=p.mesh;
    m=splitCompoundMesh(p);
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
    if typeof(p.compoundData).parameters[1]!=:nomethod
        celltype = VTKCellTypes.VTK_POLYGON;
    end
    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nf=m.topology.size[3];
    compound_nf=compound_m.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end
    fComp=Array{Array{Float64},1}(undef, length(comp))
    for l in 1:length(comp)
        fComp[l]=getCompoundElementProperties(p.femType[comp[l]][1],m.meshType,p.compoundData.assembledPhi[p.femType[comp[l]][1]],mx,my,p.compoundData);
    end
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;
    vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

    sol=p.solution[tend];
    for l in 1:length(comp)
        solc=getfield(sol,comp[l]);
        if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
            nSubCells=p.compoundData.nSubCells;
            cvtk=zeros(Float64, nf)
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                #cvtk[k]=dot(fComp[l],cLoc);
                for subCell in 1:nSubCells
                    cvtk[(k-1)*nSubCells+subCell]=mean(cLoc);
                end
            end
            vtk_cell_data(vtk, cvtk, name[l])
        else
            cvtk=zeros(Float64, m.geometry.dim, nf)
            center=Array{Float64,1}(undef,2);
            subcoord=Array{Array{Float64,2},1}(undef,p.compoundData.nSubCells);
            assembledPhi=p.compoundData.assembledPhi[p.femType[comp[l]][1]]
            divphi=p.degFBoundary[p.femType[comp[l]][1]].divphi;
            sq=length(p.compoundData.quadWeights);
            nquadPhi=p.compoundData.nquadPhi[p.femType[comp[l]][1]];
            nSubCells=p.compoundData.nSubCells;
            nquadPoints=p.compoundData.nquadPoints;
            quadWeights=p.compoundData.quadWeights;
            J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
            ddJ_edge=Array{Float64,1}(undef,sq);
            jphi_edge=initJacobi((m.geometry.dim,size(fComp[l],2)),sq);
            for i in 1:nSubCells
                #fill! causes mutating all entries of subcoord when changing a single entry
                subcoord[i]=Array{Float64,2}(undef,2,p.compoundData.nVerticesSubElement);
            end
            for k in 1:nf
                cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                coord= compound_m.geometry.coordinates[:,compound_m.topology.incidence["20"][compound_m.topology.offset["20"][k]:compound_m.topology.offset["20"][k+1]-1]]
                getSubCells!(subcoord, coord, center, p.compoundData);
                assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, quadWeights, p.compoundData);
                for subCell in 1:nSubCells
                    dJ=jacobi!(J,mx,my,subcoord[subCell],m.meshType);
                    fLoc=(1/dJ)*J*fComp[l]
                    for i in 1:length(assembledPhi)
                        for subi in 1:size(fComp[l],2)
                            cvtk[:,(k-1)*nSubCells+subCell]+=assembledPhi[i][subi,subCell]*fLoc[:,subi]*cLoc[i];
                        end
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
        end
    end

    outfiles=vtk_save(vtk);
    return outfiles::Vector{String}
end

function split_unstructured_vtk(p::femProblem, t::Array{Float64,1}, comp::Array{Symbol,1}, name::Array{String,1}, filename::String)

    compound_m=p.mesh;
    m=splitCompoundMesh(p);
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
    if typeof(p.compoundData).parameters[1]!=:nomethod
        celltype = VTKCellTypes.VTK_POLYGON;
    end
    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nf=m.topology.size[3];
    compound_nf=compound_m.topology.size[3];
    for k in 1:nf
        inds=inc[off[k]:off[k+1]-1];
        push!(cells, MeshCell(celltype, inds))
    end
    fComp=Array{Array{Float64},1}(undef, length(comp))
    for l in 1:length(comp)
        fComp[l]=getCompoundElementProperties(p.femType[comp[l]][1],m.meshType,mx,my,p.compoundData);
    end
    J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
    dJ=0.0;
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    vtk_filename_noext = (@__DIR__)*"/VTK/"*filename;

    outfiles = paraview_collection(vtk_filename_noext) do pvd
        for it in 1:length(t)
            vtk = vtk_grid(@sprintf("%s_%02i", vtk_filename_noext, it), pts, cells)
            sol=p.solution[t[it]];
            for l in 1:length(comp)
                solc=getfield(sol,comp[l]);
                if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
                    nSubCells=p.compoundData.nSubCells;
                    cvtk=zeros(Float64, nf)
                    for k in 1:compound_nf
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        #cvtk[k]=dot(fComp[l],cLoc);
                        for subCell in 1:nSubCells
                            cvtk[(k-1)*nSubCells+subCell]=mean(cLoc);
                        end
                    end
                    vtk_cell_data(vtk, cvtk, name[l])
                else
                    #=
                    cvtk=zeros(Float64, m.geometry.dim, nf)
                    for k in 1:nf
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        #dJ=jacobi!(J,m,k,mx,my,coord);
                        #fLoc=(1/dJ)*J*fComp[l]
                        #cvtk[:,k]=fLoc*cLoc;
                        cvtk[:,k]=[mean(cLoc), mean(cLoc)];
                    end
                    =#
                    cvtk=zeros(Float64, m.geometry.dim, nf)
                    center=Array{Float64,1}(undef,2);
                    subcoord=Array{Array{Float64,2},1}(undef,p.compoundData.nSubCells);
                    assembledPhi=p.compoundData.assembledPhi[p.femType[comp[l]][1]]
                    divphi=p.degFBoundary[p.femType[comp[l]][1]].divphi;
                    sq=length(p.compoundData.quadWeights);
                    nquadPhi=p.compoundData.nquadPhi[p.femType[comp[l]][1]];
                    nSubCells=p.compoundData.nSubCells;
                    nquadPoints=p.compoundData.nquadPoints;
                    quadWeights=p.compoundData.quadWeights;
                    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
                    ddJ_edge=Array{Float64,1}(undef,sq);
                    jphi_edge=initJacobi((m.geometry.dim,size(fComp[l],2)),sq);
                    for i in 1:nSubCells
                        #fill! causes mutating all entries of subcoord when changing a single entry
                        subcoord[i]=Array{Float64,2}(undef,2,p.compoundData.nVerticesSubElement);
                    end
                    for k in 1:compound_nf
                        cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
                        coord= compound_m.geometry.coordinates[:,compound_m.topology.incidence["20"][compound_m.topology.offset["20"][k]:compound_m.topology.offset["20"][k+1]-1]]
                        getSubCells!(subcoord, coord, center, p.compoundData);
                        assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, quadWeights, p.compoundData);
                        for subCell in 1:nSubCells
                            dJ=jacobi!(J,mx,my,subcoord[subCell],m.meshType);
                            fLoc=(1/dJ)*J*fComp[l]
                            for i in 1:length(assembledPhi)
                                for subi in 1:size(fComp[l],2)
                                    cvtk[:,(k-1)*nSubCells+subCell]+=assembledPhi[i][subi,subCell]*fLoc[:,subi]*cLoc[i];
                                end
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
                end
            end
            vtk_save(vtk)
            collection_add_timestep(pvd, vtk, it)
        end
    end

    return outfiles::Vector{String}
end

function compound_vtk(m::mesh, degF::degF{2}, compoundData::compoundData, sol::Array{Float64,1}, femType::Symbol, filename::String, name::String="Test")

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
    if typeof(compoundData).parameters[1]!=:nomethod
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

    #=
    for k in 1:nf
        cLoc=sol[l2g(degF, k)]
        dJ=jacobi!(J,m,k,mx,my,coord);
        fLoc=(1/dJ)*J*fComp
        cvtk[:,k]=fLoc*cLoc;
    end
    =#
    cvtk=zeros(Float64, m.geometry.dim, nf)
    fCompoundLoc=zeros(Float64, m.geometry.dim, compoundData.nCompoundPhi[femType])
    center=Array{Float64,1}(undef,2);
    subcoord=Array{Array{Float64,2},1}(undef,compoundData.nSubCells);
    assembledPhi=compoundData.assembledPhi[femType]
    divphi=degF.divphi;
    sq=length(compoundData.quadWeights);
    nquadPhi=compoundData.nquadPhi[femType];
    nSubCells=compoundData.nSubCells;
    nquadPoints=compoundData.nquadPoints;
    quadWeights=compoundData.quadWeights;
    J_edge=initJacobi((m.geometry.dim,m.topology.dim),sq);
    ddJ_edge=Array{Float64,1}(undef,sq);
    jphi_edge=initJacobi((m.geometry.dim,size(fComp,2)),sq);
    for i in 1:nSubCells
        #fill! causes mutating all entries of subcoord when changing a single entry
        subcoord[i]=Array{Float64,2}(undef,2,compoundData.nVerticesSubElement);
    end
    for k in 1:nf
        cLoc=sol[l2g(degF, k)]
        coord= m.geometry.coordinates[:,m.topology.incidence["20"][m.topology.offset["20"][k]:m.topology.offset["20"][k+1]-1]]
        getSubCells!(subcoord, coord, center, compoundData);
        assemblePhi!(assembledPhi, subcoord, m, divphi, J_edge, ddJ_edge, jphi_edge, nquadPhi, nquadPoints, quadWeights, compoundData);
        fill!(fCompoundLoc,0.0)
        for subCell in 1:nSubCells
            dJ=jacobi!(J,mx,my,subcoord[subCell],m.meshType);
            fLoc=(1/dJ)*J*fComp
            for i in 1:length(assembledPhi)
                for subi in 1:size(fComp,2)
                    fCompoundLoc[:,i]+=(1/nSubCells)*(assembledPhi[i][subi,subCell]*fLoc[:,subi]);
                end
            end
        end
        cvtk[:,k]=fCompoundLoc*cLoc;
        #println(fCompoundLoc)
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
