using TriangleMesh
import TriangleMesh.create_mesh,TriangleMesh.Mesh_ptr_C

function create_mesh(poly :: Polygon_pslg;
                                info_str :: String = "Triangular mesh of polygon (PSLG)",
                                verbose :: Bool = false,
                                check_triangulation :: Bool = false,
                                voronoi :: Bool = false,
                                delaunay :: Bool = false,
                                mesh_convex_hull :: Bool = false,
                                output_edges :: Bool = true,
                                output_cell_neighbors :: Bool = true,
                                quality_meshing :: Bool = true,
                                prevent_steiner_points_boundary :: Bool = false,
                                prevent_steiner_points :: Bool = false,
                                set_max_steiner_points :: Bool = false,
                                input_area_max :: Float64 = 0.0,
                                set_angle_min :: Bool = false,
                                add_switches :: String = "")

# this is a custom method taken and adapted from the TriangleMesh pkg
# where (if chosen to be given manually) the maximum area does not get asked for (like in the package)
# but has to be given as an optional argument

    switches = "p"

    if ~verbose
        switches = switches * "Q"
    end

    if check_triangulation
        switches = switches * "C"
    end

    if mesh_convex_hull
        switches = switches * "c"
    end

    if voronoi
        switches = switches * "v"
    end

    if delaunay
        switches = switches * "D"
    end

    if output_edges
        switches = switches * "e"
    end

    if output_cell_neighbors
        switches = switches * "n"
    end

    # -------
    if prevent_steiner_points
        prevent_steiner_points_boundary = true
        switches = switches * "YY"
    else
        if prevent_steiner_points_boundary
            switches = switches * "Y"
        end
    end

    # -------

    # -------
    if set_max_steiner_points
        max_steiner_points_str = input("Maximum number of Steiner points allowed: ")
        # Check if input makes sense
        try
            number = parse(Int, max_steiner_points_str)
            if number<0
                Base.@error("Maximum number of Steiner points must be nonnegative.")
            end
        catch
            Base.@error("Maximum number of Steiner points must be a nonnegative integer.")
        end

        switches = switches * "S" * max_steiner_points_str
    end
    # -------

    # -------
    if input_area_max!=0.0
        #max_area_str = input("Maximum triangle area: ")
        max_area_str = "$input_area_max";
        # Check if input makes sense
        try
            number = parse(Float64, max_area_str)
            if number<=0
                Base.@error "Area must be positive."
            end
        catch
            Base.@error "Area must be a positive real number."
        end

        switches = switches * "a" * max_area_str
    end
    # -------

    # -------
    if set_angle_min
        quality_meshing = true
        max_angle_str = input("Set angle constraint (choose whith care): ")
        # Check if input makes sense
        try
            number = parse(Float64, max_angle_str)
            if number <=0
                Base.@error "Minimum angle must be positive."
            end
            if number >=34
                Base.@warn "Minimum angle should not be larger than 34 degrees. For a larger angle TRIANGLE might not converge."
            end
        catch
            Base.@error "Area must be a positive real number and should not be larger than 34 degrees."
        end
        switches = switches * "q" * max_angle_str
    else
        if quality_meshing
            switches = switches * "q"
        end
    end
    # -------


    # This enables to use aditional switches and should be used with care
    switches = switches * add_switches

    if occursin("z",switches)
        Base.@error("Triangle switches must not contain `z`. Zero based indexing is not allowed.")
    end

    mesh_in = Mesh_ptr_C(poly)

    mesh_buffer = Mesh_ptr_C()
    vor_buffer = Mesh_ptr_C()

    triangulate(mesh_in, mesh_buffer, vor_buffer, switches)
    mesh = TriMesh(mesh_buffer, vor_buffer, info_str)

    return mesh
end
