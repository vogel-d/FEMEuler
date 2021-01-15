function splitCompoundMeshKites(p::femProblem)
    compoundMesh=p.mesh;
    size=compoundMesh.topology.size;
    nSubCells=p.data.compoundData.nSubCells;
    mt=compoundMesh.meshType;
    n=copy(compoundMesh.topology.n);
    condEW=compoundMesh.boundaryConditionEW;
    condTB=compoundMesh.boundaryConditionTB;
    bE=spzeros(Int, size[2]*2+size[3]*nSubCells);
    bV=spzeros(Int, size[1]+size[2]+size[3]);
    compoundbE=compoundMesh.boundaryEdges;
    compoundbV=compoundMesh.boundaryVertices;
    xl=compoundMesh.geometry.l[1]; yl=compoundMesh.geometry.l[2];
    xr=compoundMesh.geometry.r[1]; yr=compoundMesh.geometry.r[2];
    inc20=compoundMesh.topology.incidence["20"];
    inc10=compoundMesh.topology.incidence["10"];
    inc01=compoundMesh.topology.incidence["01"];
    off20=compoundMesh.topology.offset["20"];
    off10=compoundMesh.topology.offset["10"];
    off01=compoundMesh.topology.offset["01"];
    String(typeof(p.data.compoundData).parameters[1])[(end-4):end]!="Kites" && error("splitting available for kites only")


    coordx=copy(compoundMesh.geometry.coordinates[1,:]);
    coordy=copy(compoundMesh.geometry.coordinates[2,:]);
    for edge in 1:size[2]
        coord_edge=compoundMesh.geometry.coordinates[:,inc10[off10[edge]:off10[edge+1]-1]]
        push!(coordx,0.5*(coord_edge[1,1]+coord_edge[1,2]))
        push!(coordy,0.5*(coord_edge[2,1]+coord_edge[2,2]))
    end
    for face in 1:size[3]
        coord_face=compoundMesh.geometry.coordinates[:,inc20[off20[face]:off20[face+1]-1]]
        push!(coordx,mean(coord_face[1,:]))
        push!(coordy,mean(coord_face[2,:]))
    end

    coord=Array{Float64,2}(undef,2,size[1]+size[2]+size[3]);
    coord[1,:]=coordx;
    coord[2,:]=coordy;

    inc20_new=Int64[];
    off20_new=sort(collect(range(1,step=mt,length=1+size[3]*nSubCells)));
    inner_edges=Int64[];
    for row in 1:n[2]
        for col in 1:n[1]
            Cell=(row-1)*n[1]+col;
            vertices=inc20[off20[Cell]:off20[Cell+1]-1];
            for subCell in 1:3
                #detect compound vertices that define the specific subCell
                vertex1=vertices[mod(subCell+nSubCells-2,nSubCells)+1]
                vertex2=vertices[subCell]
                vertex3=vertices[mod(subCell,nSubCells)+1]
                #get the edge number (not ordered in inc) to detect the new vertex indices
                edge1=intersect(inc01[off01[vertex1]:off01[vertex1+1]-1],inc01[off01[vertex2]:off01[vertex2+1]-1])[1];
                edge2=intersect(inc01[off01[vertex2]:off01[vertex2+1]-1],inc01[off01[vertex3]:off01[vertex3+1]-1])[1];

                #push kite-subcell
                #append!(inc20_new,[size[1]+edge1,vertex2,size[1]+edge2,size[1]+size[2]+Cell]);
                append!(inc20_new,[vertex2,size[1]+edge2,size[1]+size[2]+Cell,size[1]+edge1]);

                #push inner edge
                append!(inner_edges,[size[1]+edge2,size[1]+size[2]+Cell])
            end
        end
        for col in 1:n[1]
            Cell=(row-1)*n[1]+col;
            vertices=inc20[off20[Cell]:off20[Cell+1]-1];
            for subCell in 6:-1:4
                #detect compound vertices that define the specific subCell
                vertex1=vertices[mod(subCell+nSubCells-2,nSubCells)+1]
                vertex2=vertices[subCell]
                vertex3=vertices[mod(subCell,nSubCells)+1]
                #get the edge number (not ordered in inc) to detect the new vertex indices
                edge1=intersect(inc01[off01[vertex1]:off01[vertex1+1]-1],inc01[off01[vertex2]:off01[vertex2+1]-1])[1];
                edge2=intersect(inc01[off01[vertex2]:off01[vertex2+1]-1],inc01[off01[vertex3]:off01[vertex3+1]-1])[1];

                #push kite-subcell
                #append!(inc20_new,[size[1]+edge1,vertex2,size[1]+edge2,size[1]+size[2]+Cell]);
                append!(inc20_new,[vertex2,size[1]+edge2,size[1]+size[2]+Cell,size[1]+edge1]);

                #push inner edge
                append!(inner_edges,[size[1]+edge2,size[1]+size[2]+Cell])
            end
        end
    end

    inc10_new=Int64[];
    off10_new=sort(collect(range(1,step=2,length=1+size[2]*2+size[3]*nSubCells)));
    for edge in 1:size[2]
        vertices=inc10[off10[edge]:off10[edge+1]-1];
        #sort! for continuity in ordering at boundary
        sort!(vertices);
        append!(inc10_new,[vertices[1],size[1]+edge]);
        append!(inc10_new,[size[1]+edge,vertices[2]]);
        if compoundbE[edge]==1.0
            bE[edge*2-1]    =1.0;
            bE[edge*2]      =1.0;
            bV[size[1]+edge]=1.0;
        elseif compoundbE[edge]<0.0
            bE[edge*2-1]=-(abs(compoundbE[edge])*2-1);
            bE[edge*2]  =-(abs(compoundbE[edge])*2);
            bV[size[1]+edge]=-(size[1]+abs(compoundbE[edge]))
        end
    end

    bV[1:size[1]]=compoundbV;
    append!(inc10_new,inner_edges);

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>inc20_new,"10"=>inc10_new);
    off=Dict("20"=>off20_new,"10"=>off10_new);

    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n[1]=n[1]*3;
    n[2]=n[2]*2;
    l=Float64[xl,yl];
    r=Float64[xr,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV,condEW,condTB);

    return m
end



function splitCompoundMeshTris(p::femProblem)
    compoundMesh=p.mesh;
    size=compoundMesh.topology.size;
    nSubCells=p.data.compoundData.nSubCells;
    mt=compoundMesh.meshType;
    condEW=compoundMesh.boundaryConditionEW;
    condTB=compoundMesh.boundaryConditionTB;
    bE=spzeros(Int, size[2]*2+size[3]*nSubCells);
    bV=spzeros(Int, size[1]+size[2]+size[3]);
    compoundbE=compoundMesh.boundaryEdges;
    compoundbV=compoundMesh.boundaryVertices;
    nCompoundVertices=diff(compoundMesh.topology.offset["20"][1:2])[1]
    xl=compoundMesh.geometry.l[1]; yl=compoundMesh.geometry.l[2];
    xr=compoundMesh.geometry.r[1]; yr=compoundMesh.geometry.r[2];
    inc20=compoundMesh.topology.incidence["20"];
    inc10=compoundMesh.topology.incidence["10"];
    inc01=compoundMesh.topology.incidence["01"];
    off20=compoundMesh.topology.offset["20"];
    off10=compoundMesh.topology.offset["10"];
    off01=compoundMesh.topology.offset["01"];

    coordx=copy(compoundMesh.geometry.coordinates[1,:]);
    coordy=copy(compoundMesh.geometry.coordinates[2,:]);
    for edge in 1:size[2]
        coord_edge=compoundMesh.geometry.coordinates[:,inc10[off10[edge]:off10[edge+1]-1]]
        push!(coordx,0.5*(coord_edge[1,1]+coord_edge[1,2]))
        push!(coordy,0.5*(coord_edge[2,1]+coord_edge[2,2]))
    end
    for face in 1:size[3]
        coord_face=compoundMesh.geometry.coordinates[:,inc20[off20[face]:off20[face+1]-1]]
        push!(coordx,mean(coord_face[1,:]))
        push!(coordy,mean(coord_face[2,:]))
    end

    coord=Array{Float64,2}(undef,2,size[1]+size[2]+size[3]);
    coord[1,:]=coordx;
    coord[2,:]=coordy;

    inc20_new=Int64[];
    off20_new=sort(collect(range(1,step=mt,length=1+size[3]*nSubCells)));
    inner_edges=Int64[];
    for Cell in 1:size[3]
        vertices=inc20[off20[Cell]:off20[Cell+1]-1];
        for subCell in 1:nCompoundVertices
            #detect compound vertices that define the specific subCell
            vertex1=vertices[mod(subCell+nCompoundVertices-2,nCompoundVertices)+1]
            vertex2=vertices[subCell]
            vertex3=vertices[mod(subCell,nCompoundVertices)+1]
            #get the edge number (not ordered in inc) to detect the new vertex indices
            edge1=intersect(inc01[off01[vertex1]:off01[vertex1+1]-1],inc01[off01[vertex2]:off01[vertex2+1]-1])[1];
            edge2=intersect(inc01[off01[vertex2]:off01[vertex2+1]-1],inc01[off01[vertex3]:off01[vertex3+1]-1])[1];

            #push kite-subcell
            append!(inc20_new,[vertex2,size[1]+edge1,size[1]+size[2]+Cell])
            append!(inc20_new,[vertex2,size[1]+edge2,size[1]+size[2]+Cell])

            #push inner edge
            append!(inner_edges,[size[1]+edge2,size[1]+size[2]+Cell])
            append!(inner_edges,[vertex2,size[1]+size[2]+Cell])
        end
    end

    inc10_new=Int64[];
    off10_new=sort(collect(range(1,step=2,length=1+size[2]*2+size[3]*nSubCells)));
    for edge in 1:size[2]
        vertices=inc10[off10[edge]:off10[edge+1]-1];
        #sort! for continuity in ordering at boundary
        sort!(vertices);
        append!(inc10_new,[vertices[1],size[1]+edge]);
        append!(inc10_new,[size[1]+edge,vertices[2]]);
        if compoundbE[edge]==1.0
            bE[edge*2-1]    =1.0;
            bE[edge*2]      =1.0;
            bV[size[1]+edge]=1.0;
        elseif compoundbE[edge]<0.0
            bE[edge*2-1]=-(abs(compoundbE[edge])*2-1);
            bE[edge*2]  =-(abs(compoundbE[edge])*2);
            bV[size[1]+edge]=-(size[1]+abs(compoundbE[edge]))
        end
    end

    bV[1:size[1]]=compoundbV;
    append!(inc10_new,inner_edges);

    #Initialisieren der Inzidenz mit den Einträgen "20" und "10"
    inc=Dict("20"=>inc20_new,"10"=>inc10_new);
    off=Dict("20"=>off20_new,"10"=>off10_new);

    #Initialisieren der Topologie, Geometrie und damit des Meshes
    n=compoundMesh.topology.n;
    l=Float64[xl,yl];
    r=Float64[xr,yr];
    mT=meshTopology(inc,off,n);
    mG=meshGeometry(coord,l,r);
    m=mesh(mT,mG, bE, bV,condEW,condTB);

    return m
end
