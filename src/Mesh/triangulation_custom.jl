include("create_mesh_custom.jl");

function triangles(x_range::Array{Float64,1},y_range::Array{Float64,1}, rows::Int64, condEW::Symbol, condTB::Symbol)

    ###set domain

    cols=Int(ceil((diff(x_range)[1]/diff(y_range)[1])*rows));
    #stepwidth
    dy=diff(y_range)[1]/rows;
    dx=diff(x_range)[1]/cols;

    #domain vertices (with boundary vertices)
    domain_nodes=Array{Float64,2}(undef,2*cols+2*rows,2);
    vertex_id=1
    for i in 0:cols-1
        domain_nodes[vertex_id,:]=[x_range[1]+dx*i,y_range[1]];
        vertex_id+=1;
    end
    for i in 0:rows-1
        domain_nodes[vertex_id,:]=[x_range[2],y_range[1]+dy*i];
        vertex_id+=1;
    end
    for i in cols:-1:1
        domain_nodes[vertex_id,:]=[x_range[1]+dx*i,y_range[2]];
        vertex_id+=1;
    end
    for i in rows:-1:1
        domain_nodes[vertex_id,:]=[x_range[1],y_range[1]+dy*i];
        vertex_id+=1;
    end

    #domain boundary
    domain_boundary=Array{Int64,2}(undef,2*rows+2*cols,2);
    for i in 1:2*rows+2*cols
        domain_boundary[i,1]=i;
        domain_boundary[i,2]=i+1;
    end
    domain_boundary[end,2]=1;

    #classifying (all same)
    #node_marker=ones(Int,4,1);
    #boundary_marker=2*ones(Int,size(domain_boundary,1));

    domain=Polygon_pslg(2*cols+2*rows,1,1,2*rows+2*cols,0);

    set_polygon_point!(domain,domain_nodes);
    set_polygon_segment!(domain,domain_boundary);
    #set_polygon_segment_marker!(domain, boundary_marker)

    #determine resolution
    #pick max_area by choosing a little more than the area
    #of an equilateral triangle of edge length dy
    max_area=1.2*((sqrt(3)/4)*dy*dy);

    #mesh_raw = create_mesh(domain, info_str="my mesh", voronoi=true, delaunay=true, set_area_max=true, prevent_steiner_points_boundary=true);
    mesh_raw = create_mesh(domain, info_str="my mesh", voronoi=true, delaunay=true, input_area_max=max_area, prevent_steiner_points_boundary=true);

    inc20=Int64[];
    for i in 1:mesh_raw.n_cell
        #sort for correct local numbering
        append!(inc20,sort(mesh_raw.cell[:,i]));
    end
    inc10=Int64[];
    for i in 1:mesh_raw.n_edge
        #sort for correct local numbering, maybe not necessary
        append!(inc10,sort(mesh_raw.edge[:,i]));
    end
    inc=Dict("20"=>inc20,"10"=>inc10);

    off20=sort(collect(1:3:(mesh_raw.n_cell*3+1)));
    off10=sort(collect(1:2:(mesh_raw.n_edge*2+1)));
    off=Dict("20"=>off20,"10"=>off10);

    #default irrelevant n
    n=[0,0];

    topology=meshTopology(inc,off,n);
    geometry=meshGeometry(mesh_raw.point,[x_range[1],y_range[1]],[x_range[2],y_range[2]]);


    bV=spzeros(Int, topology.size[1]);
    bE=spzeros(Int, topology.size[2]);

    #set boundaryVertices
    corner_vertices=[1,cols+1,cols+rows+1,2*cols+rows+1];
    left_side_vertices=  sort!(collect((2*cols+rows+2):(2*rows+2*cols)));
    right_side_vertices= sort!(collect((cols+2):(cols+rows)),rev=true);
    upper_side_vertices= sort!(collect((cols+rows+2):(2*cols+rows)),rev=true);
    down_side_vertices=  sort!(collect(2:cols));

    if condEW==:periodic
        if condTB==:periodic
            bV[corner_vertices[2:4]].=-1;
            bV[down_side_vertices]=-upper_side_vertices;
            bV[left_side_vertices]=-right_side_vertices;
        elseif condTB==:constant
            bV[corner_vertices].=1;
            bV[down_side_vertices].=1;
            bV[upper_side_vertices].=1;
            bV[left_side_vertices]=-right_side_vertices;
        end
    elseif condEW==:constant
        if condTB==:periodic
            bV[corner_vertices].=1;
            bV[left_side_vertices].=1;
            bV[right_side_vertices].=1;
            bV[down_side_vertices]=-upper_side_vertices;
        elseif condTB==:constant
            bV[corner_vertices].=1;
            bV[right_side_vertices].=1;
            bV[left_side_vertices].=1;
            bV[down_side_vertices].=1;
            bV[upper_side_vertices].=1;
        end
    end

    #set boundaryEdges
    left_side_vertices_big=  vcat(left_side_vertices,corner_vertices[[4,1]]);
    right_side_vertices_big= vcat(right_side_vertices,corner_vertices[[2,3]]);
    upper_side_vertices_big= vcat(upper_side_vertices,corner_vertices[[3,4]]);
    down_side_vertices_big=  vcat(down_side_vertices,corner_vertices[[1,2]]);
    for edge_id in 1:size(mesh_raw.edge,2)
        point_1=mesh_raw.edge[1,edge_id];
        point_2=mesh_raw.edge[2,edge_id];
        #edge at down side?
        if in(point_1,down_side_vertices_big) && in(point_2,down_side_vertices_big)
            if condTB==:periodic
                #detect opposite vertices
                opposite_vertices=[point_1+2*(cols+1-point_1)+rows;
                                   point_2+2*(cols+1-point_2)+rows];

                #find opposite edge by searching for opposite vertices
                case_1=findall(opposite_vertices,mesh_raw.edge);
                case_2=findall(reverse!(opposite_vertices),mesh_raw.edge);
                opposite_edge=append!(case_1,case_2)[1];
                bE[edge_id]=-opposite_edge;
            elseif condTB==:constant
                bE[edge_id]=1;
            end
        #edge at right side?
        elseif in(point_1,right_side_vertices_big) && in(point_2,right_side_vertices_big)
            if condEW==:periodic
                #nothing needed for right side
                continue;
            elseif condEW==:constant
                bE[edge_id]=1;
            end
        #edge at upper side?
        elseif in(point_1,upper_side_vertices_big) && in(point_2,upper_side_vertices_big)
            if condTB==:periodic
                #nothing needed for upper side
                continue;
            elseif condTB==:constant
                bE[edge_id]=1;
            end
        #edge at left side?
        elseif in(point_1,left_side_vertices_big) && in(point_2,left_side_vertices_big)
            if condEW==:periodic
                #detect opposite vertices
                    #handle inconsistent numbering at left side
                    if point_1==1 point_1=2*cols+2*rows+1 end;
                    if point_2==1 point_2=2*cols+2*rows+1 end;
                    #renumber, counted from  bottom left corner
                    point_1=2*cols+2*rows+2-point_1;
                    point_2=2*cols+2*rows+2-point_2;
                opposite_vertices=[cols+point_1;
                                   cols+point_2]

                #find opposite edge by searching for opposite vertices
                case_1=findall(opposite_vertices,mesh_raw.edge);
                case_2=findall(reverse!(opposite_vertices),mesh_raw.edge);
                opposite_edge=append!(case_1,case_2)[1];
                bE[edge_id]=-opposite_edge;
            elseif condEW==:constant
                bE[edge_id]=1;
            end
        end
    end

    m=mesh(topology,geometry,bE,bV,condEW,condTB);

    return m;
end
