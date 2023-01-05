module WWUMeshes

    # Exported Functions
    export read_obj, write_obj
    export gen_mesh, est_normals, cmp_mesh, compress_obj
    export cube_mesh, cylinder_mesh, sphere_mesh, estimate_normals
    export OBJTriangle, OBJMesh

    # Usings & Settings
    using FileIO
    using LinearAlgebra

    push!(LOAD_PATH, pwd())

    include("GfxBase.jl")
    using .GfxBase


    # OBJTriangle - struct that represents a single triangle in a mesh.
    mutable struct OBJTriangle
        positions::Array{Int, 1} # vertex position indices
        uvs::Array{Int, 1} # vertex texture coordinate indices
        normals::Array{Int, 1} # normal vector indices
    end

    #OBJMesh - struct that represents an indexed triangle mesh for reading from or writing to OBJ format. 
    mutable struct OBJMesh
        positions::Array{Vec3, 1} # all vertex positions
        uvs::Array{Vec2, 1} # all texture coordinates
        normals::Array{Vec3, 1} # all vertex normals
        triangles::Array{OBJTriangle, 1} # the OBJTriangles belonging to the mesh
    end

   #= read_obj(obj_filename)
    = 
    = Read a mesh in OBJ format from file obj_filename.
    =#
    function read_obj(obj_filename)
        m = OBJMesh([], [], [], []) # create a mesh
        open(obj_filename) do f
            for (line_number, line) in enumerate(eachline(f))
                if line == "" || line[1] == "#"
                    continue # skip comments
                end
                # Read the line and add its contents to the correct field of m:
                tokens = split(strip(line))
                if tokens[1] == "v" # vertex
                    push!(m.positions, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
                elseif tokens[1] == "vt" # vertex texture
                    push!(m.uvs, Vec2([parse(Float64, x) for x in tokens[2:end]]...))
                elseif tokens[1] == "vn" # vertex normal
                    push!(m.normals, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
                elseif tokens[1] == "f"
                    # create a OBJTriangle face:
                    points = []
                    uvs = []
                    normals = []
                    # handle faces with no texture and/or normals
                    for corner in tokens[2:end]
                        indices = split(corner, '/')
                        if length(indices) == 3 # all 3 present, third is a normal
                            push!(normals, parse(Int, indices[3]))
                        end
                        if length(indices) >= 2 && indices[2] != ""
                            # if there are 2 or more and the second isn't blank, it's a texture
                            push!(uvs, parse(Int, indices[2]))
                        end
                        if length(indices) >= 1 # first value is the position
                            push!(points, parse(Int, indices[1]))
                        else # unless it has none, in which case it's not valid
                            error("in line $line_number: face vertex $corner could not be parsed")
                        end
                    end
                    # create the triangle and add it to the triangles array
                    push!(m.triangles, OBJTriangle(points, uvs, normals))
                end
            end
        end
        return m
    end

    """ write_obj(obj_filename)
    Write the given mesh in OBJ format to file obj_filename."""
    function write_obj(obj_filename, mesh::OBJMesh)
        open(obj_filename, "w") do f
            # write all positions:
            for v in mesh.positions
                write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
            end

            # write all texture coords:
            for v in mesh.uvs
                write(f, "vt $(v[1]) $(v[2])\n")
            end
            # write all normals:
            for v in mesh.normals
                write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
            end

            # write all triangles:
            for tri in mesh.triangles
                write(f, "f $(tri_vertex_str(tri))\n")
            end

        end

    end

   #= tri_vertex_str(triangle)
    = 
    = Return a string with the indices of applicable positions, texture coordinates,
    = and normals for a given triangle according to the OBJ specification.
    = In particular, if p, u, and n are position, vertex and normal, each corner
    = of the triangle is represented as one of the following:
    =     p       (position only)
    =     p/u     (position and texture)
    =     p//n    (position and normal)
    =     p/u/n   (position, texture, and normal) 
    =#
    function tri_vertex_str(triangle::OBJTriangle)
        # determine whether textures and normals are present:
        write_uv = length(triangle.uvs) == length(triangle.positions)
        write_normals = length(triangle.normals) == length(triangle.positions)
        corners = []
        for i = 1:3
            output = "$(triangle.positions[i])"
            if write_uv && !write_normals
                output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
            elseif !write_uv && write_normals
                output = output * "//$(triangle.normals[i])"
            elseif write_uv && write_normals
                output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
            end
            push!(corners, output)
        end
        join(corners, " ")
    end


    #=
    = gen_mesh(outfile, geom, divisionsU, divisionsV)
    = 
    = Generate a mesh and save the result in a file with name outfile.
    = geom may be "cube", "cylinder", "sphere", or "torus".
    = Cube requires nothing; cylinder requires divisionsU; sphere requires divisionsU and divisionsV; torus requires divisions, divisionsV, and minorRadius.
    =#
    function gen_mesh(outfile, geom, divisionsU=32, divisionsV=16, minorRadius=0.2)
        if geom == "cube"
            mesh = cube_mesh()
        elseif geom == "cylinder"
            mesh = cylinder_mesh(divisionsU)
        elseif geom == "sphere"
            mesh = sphere_mesh(divisionsU, divisionsV)
        elseif geom == "torus"
            mesh = torus_mesh(divisionsU, divisionsV, minorRadius);
        end
        write_obj(outfile, mesh)
    end


   #=
    = est_normals(outfile, infile)
    =
    = Estimate normals of the mesh stored in infile, saving the result in outfile.
    =#
    function est_normals(outfile, infile)
        input_mesh = read_obj(infile)
        mesh = estimate_normals(input_mesh)
        write_obj(outfile, mesh)
    end

   #=
    = cmp_mesh(infile1, infile2, verbose, epsilon)
    =
    = Determine if 2 meshes m1, m2 are equivalent under the following definition of equivalency:
        =   1. The two meshes have the same number of faces.
        =   2. For each face in m1, a face exists in m2 with the same vertex positions, texture coordinates, and normals (if applicable)
        =   3. For each face in m2, a face exists in m1 with the same vertex positions, texture coordinates, and normals (if applicable)
    =#
    function cmp_mesh(infile1, infile2, verbose, epsilon)
        input_mesh1 = read_obj(infile1);
        input_mesh2 = read_obj(infile2);
        println(mesh_compare(input_mesh1, input_mesh2, verbose, epsilon));
    end

   #=
    = cube_mesh()
    = 
    = Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
    = axis-aligned. 
    =#
    function cube_mesh()
        positions = []
        uvs = []
        normals = []
        triangles = []
        # key to comments:
        # L/R = x = right/left
        # B/T = y = top/bottom
        # C/F = z = close/far
        push!(positions, Vec3( 1, -1, -1)) # 1 RBC
        push!(positions, Vec3( 1, -1,  1)) # 2 RBF
        push!(positions, Vec3(-1, -1,  1)) # 3 LBF
        push!(positions, Vec3(-1, -1, -1)) # 4 LBC
        push!(positions, Vec3( 1,  1, -1)) # 5 RTC
        push!(positions, Vec3( 1,  1,  1)) # 6 RTF
        push!(positions, Vec3(-1,  1,  1)) # 7 LTF
        push!(positions, Vec3(-1,  1, -1)) # 8 LTC

        # texture coordinates:
        push!(uvs, Vec2(1, 1)) # TR
        push!(uvs, Vec2(0, 1)) # TL
        push!(uvs, Vec2(0, 0)) # BL
        push!(uvs, Vec2(1, 0)) # BR

        # normals:
        push!(normals, Vec3( 1, 0, 0)) # R
        push!(normals, Vec3(-1, 0, 0)) # L
        push!(normals, Vec3( 0, 1, 0)) # U
        push!(normals, Vec3( 0,-1, 0)) # D
        push!(normals, Vec3( 0, 0, 1)) # C
        push!(normals, Vec3( 0, 0,-1)) # F

        # 8 faces, 2 triangles each
        push!(triangles, OBJTriangle([1,2,3], [1,2,3], [4,4,4])) # bottom face 1
        push!(triangles, OBJTriangle([1,3,4], [1,3,4], [4,4,4])) # bottom face 2
        push!(triangles, OBJTriangle([1,5,6], [4,1,2], [1,1,1])) # right face 1
        push!(triangles, OBJTriangle([1,6,2], [4,2,3], [1,1,1])) # right face 2
        push!(triangles, OBJTriangle([2,6,7], [4,1,2], [5,5,5])) # far face 1
        push!(triangles, OBJTriangle([2,7,3], [4,2,3], [5,5,5])) # far face 2
        push!(triangles, OBJTriangle([3,7,8], [2,3,4], [2,2,2])) # left face 1
        push!(triangles, OBJTriangle([3,8,4], [2,4,1], [2,2,2])) # left face 2
        push!(triangles, OBJTriangle([4,8,5], [2,3,4], [6,6,6])) # far face 1
        push!(triangles, OBJTriangle([4,5,1], [2,4,1], [6,6,6])) # far face 2
        push!(triangles, OBJTriangle([5,8,7], [1,2,3], [3,3,3])) # top face 1
        push!(triangles, OBJTriangle([5,7,6], [1,3,4], [3,3,3])) # top face 2

        # julia automatically returns the last value in the function:
        OBJMesh(positions, uvs, normals, triangles)

    end


   #= cylinder_mesh(n)
    = 
    = Return a new OBJMesh object approximation of a cylinder with radius 1 and
    = height 2, centered at the origin. The logitudinal axis is aligned with y, and
    = it is tesselated with n divisions arranged radially around the outer surface.
    = The ends of the cylinder are disc-shaped caps parallel to the xz plane. See the
    = assignment writeup for a diagram and details.
    =#
    function cylinder_mesh(divisionsU)
        # TODO
        positions = []
        uvs = []
        normals = []
        triangles = []
    
        # known values
        v0 = Vec3(0, 1, 0)
        top_cap_norms = [1, 1, 1]
        bot_cap_norms = [2, 2, 2]
    
        push!(positions, v0) # center of top cap in x,y,z
        push!(positions, -v0) # center of bottom cap in x,y,z
        
        push!(normals, v0)  # upward pointing normal
        push!(normals, -v0) # downward pointing normal
    
        push!(uvs, Vec2(0.75, 0.75)) # center of top cap in u, v
        push!(uvs, Vec2(0.25, 0.75)) # center of bottom cap in u, v
        push!(uvs, Vec2(1, 0.5)) # upper end of the side label in u, v
        push!(uvs, Vec2(1, 0)) # lower end of side label in u, v
    
        # generate range of thetas for top/bottom caps
        thetas = range(0, 2 * pi - 2 * pi / divisionsU, length=divisionsU)
    
        for theta in thetas
            # get the x, z coords for vertices (y is known to be +-1)
            x = -sin(theta)
            z = -cos(theta)
            
            # calculate the texture coords
            # get texture u-coord for shell from theta (v is known to be either 0 or 0.5 for shell)
            u_shell = theta / (2 * pi)
            # get top cap texture coords
            u_top = 0.25 * cos(theta + pi/2) + 0.75
            v_top = 0.25 * sin(theta + pi/2) + 0.75
            # get bottom cap texture coords
            u_bottom = 0.25 * cos(theta - pi/2) + 0.25
            v_bottom = 0.25 * sin(theta - pi/2) + 0.75
    
            # push top-cap vertices to list
            push!(positions, Vec3(x, 1, z))
            # push bottom-cap vertices to list
            push!(positions, Vec3(x, -1, z))
            
            # push new normal to list
            push!(normals, Vec3(x, 0, z))
            
            # push side-label texture coords to list
            push!(uvs, Vec2(u_shell, 0.5))
            push!(uvs, Vec2(u_shell, 0))
            push!(uvs, Vec2(u_top, v_top))
            push!(uvs, Vec2(u_bottom, v_bottom))
    
            # push top triangles to list
            # top-cap triangle indices: 1, length(positions)-1, length(positions) + 1
            top_cap_vertices = [1, length(positions) - 1, length(positions) + 1]
            top_cap_textures = [1, length(uvs) - 1, length(uvs) + 3]
            push!(triangles, OBJTriangle(top_cap_vertices, top_cap_textures, top_cap_norms))
    
            # push bottom triangles to list
            # bottom-cap triangle vertices: 2, length(positions) + 2, length(positions)
            bot_cap_vertices = [2, length(positions) + 2, length(positions)]
            bot_cap_textures = [2, length(uvs) + 4, length(uvs)]
            push!(triangles, OBJTriangle(bot_cap_vertices, bot_cap_textures, bot_cap_norms))
    
            # push upward pointing triangles to list
            # upward pointing side-triangle vertices: length + 1, length, length + 2
            upward_vertices = [length(positions) + 1, length(positions), length(positions) + 2]
            upward_textures = [length(uvs)+1, length(uvs) - 2, length(uvs) + 2]
            upward_norms = [length(normals) + 1, length(normals), length(normals) + 1]
            push!(triangles, OBJTriangle(upward_vertices, upward_textures, upward_norms))
    
            # push downward pointing triangles to list
            # downward pointing side-triangle vertices: length, length + 1, length - 1
            downward_vertices = [length(positions), length(positions) + 1, length(positions) - 1]
            downward_textures = [length(uvs)-2, length(uvs) + 1, length(uvs) - 3]
            downward_norms = [length(normals), length(normals) + 1, length(normals)]
            push!(triangles, OBJTriangle(downward_vertices, downward_textures, downward_norms))
        end
    
        # correct triangle indices
    
        # correct top-cap triangle indices
        triangles[length(triangles)-3].positions[3] = 3
        triangles[length(triangles)-3].uvs[3] = 7
        # correct bottom-cap triangle indices
        triangles[length(triangles)-2].positions[2] = 4
        triangles[length(triangles)-2].uvs[2] = 8
    
        # correct upward facing triangle indices
        triangles[length(triangles)-1].positions[1] = 3
        triangles[length(triangles)-1].positions[3] = 4
        triangles[length(triangles)-1].normals[1] = 3
        triangles[length(triangles)-1].normals[3] = 3
        triangles[length(triangles)-1].uvs[1] = 3
        triangles[length(triangles)-1].uvs[3] = 4
    
        #correct downward facing triangle indices
        triangles[length(triangles)].positions[2] = 3
        triangles[length(triangles)].normals[2] = 3
        triangles[length(triangles)].uvs[2] = 3
        
        OBJMesh(positions, uvs, normals, triangles)
    end
   #= 
    = sphere_mesh(n, m)
    =
    = Creates a Latitude-Longitude-tesselated approximation of a sphere with radius 1
    = centered at the origin. There are n divisions around the equator and m
    = divisions from pole to pole along each line of longitude. The North pole is at
    = (0,1,0), and the South pole is at (0,-1,0), and points on the Greenwich meridian are
    = in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
    = with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
    = coordinate varies with latitude with v=0 at the South pole and v=1 at the North
    = pole. Normals are normal to the ideal sphere surface.
    =#
    function sphere_mesh(n, m)
        # Create mesh data arrays + poles
        textures = [];
        triangles = [];
        north = Vec3(0, 1, 0);
        south = Vec3(0, -1, 0);
        positions = [north, south];

        # Unwrapped index lookup tables. Vertices are stored as they will appear in the UV texture (with the poles as an exception).
        vert_indices = zeros(Int32, n+1, m-1);
        text_indices = zeros(Int32, n+1, m-1);
        

        # Push UVs for the triangles at the north/south poles.
        for lo = 0:(n-1)
            push!(textures, Vec2((lo)/n, 1));
        end
    

        for lo = 0:(n-1)
            push!(textures, Vec2((lo+1)/n, 0));
        end

        # Generate vertices on the sphere
        for lo = 0:n

            rad_lo = ((lo/n)*2*pi)-pi;
            u = (lo/n);

            for la = 1:(m-1)
                # Calculate the index this vertice/texcoord will have in the 'flat' array
                vert_index = (lo)*(m-1) + la + 2;
                text_index = (lo)*(m-1) + la + 2*n;

                rad_la = ((la/m)*pi) - (pi/2);
                v = la / m;

                # Generate the point using its latitude & longitude
                gen_p = Vec3(cos(rad_la)*sin(rad_lo), sin(rad_la), cos(rad_la) * cos(rad_lo));

                # Push vertex and texcoord
                push!(positions, gen_p);
                push!(textures, Vec2(u, v));

                # Store the indices of the position and texcoord for later lookup
                vert_indices[lo+1, la] = vert_index;
                text_indices[lo+1, la] = text_index;
            end
        end

        # Generate triangles on the sphere
        for x = 1:n
            for y = 1:(m-2)
                # Calculate x+1 module n to assure triangles can find vertices across the seam
                xp1 = (x + 1 > n ? x + 1 - n : x + 1);

                # Use the lookup tables to find the true vertex/texcoord coordinates
                triangle = OBJTriangle(
                    [vert_indices[x,y+1], vert_indices[x, y], vert_indices[xp1,y]], 
                    [text_indices[x,y+1], text_indices[x, y], text_indices[x+1,y]],
                    [vert_indices[x,y+1], vert_indices[x, y], vert_indices[xp1,y]]
                );
                triangle_alt = OBJTriangle(
                    [vert_indices[x,y+1], vert_indices[xp1, y], vert_indices[xp1,y+1]], 
                    [text_indices[x,y+1], text_indices[x+1, y], text_indices[x+1,y+1]], 
                    [vert_indices[x,y+1], vert_indices[xp1, y], vert_indices[xp1,y+1]]
                );

                # Push the triangles
                push!(triangles, triangle, triangle_alt);
            end
        end

        # Generate the triangles directly around the poles
        for x = 1:n
            xp1 = (x + 1 > n ? (x + 1 - n) : (x + 1));

            # Indexes of the first and last rings (of vertices)
            y_n = m-1;
            y_s = 1;

            # Create the triangles
            triangle = OBJTriangle(
                [vert_indices[x,y_n], vert_indices[xp1, y_n], 1], 
                [text_indices[x,y_n], text_indices[x+1, y_n], x], 
                [vert_indices[x,y_n], vert_indices[xp1, y_n], 1]
            );
            triangle_alt = OBJTriangle(
                [vert_indices[xp1,y_s], vert_indices[x, y_s], 2], 
                [text_indices[x+1,y_s], text_indices[x, y_s], x + n], 
                [vert_indices[xp1,y_s], vert_indices[x, y_s], 2]
            );

            # Push the triangles
            push!(triangles, triangle, triangle_alt);
        end

        # Create and return the mesh
        return OBJMesh(positions, textures, copy(positions), triangles);
    end


   #=
    = torus_mesh(n, m, r)
    =
    = Creates a tesselated approximation of a torus. Has m divisions around its major radius (R = 1) and m divisions around its minor radius (r).
    = The UV seam is located on the negative z half of the xz plane.
    =#
    function torus_mesh(n, m, r) :: OBJMesh
        # Create mesh data arrays
        textures = [];
        triangles = [];
        positions = [];
        normals = [];

        # Unwrapped index lookup tables. Vertices are stored as they will appear in the UV texture.
        vert_indices = zeros(Int32, n, m);
        text_indices = zeros(Int32, n+1, m+1);
        
        # Generate points on the torus
        for lo = 0:(n-1)

            rad_lo = ((lo/n)*2*pi)-pi;
            u = lo/n;

            # Calculate a point inside the torus at the current longitude
            inner_point = Vec3(sin(rad_lo), 0, cos(rad_lo));

            for la = 0:(m-1)
                # Calculate the index this vertice/texcoord will have in the 'flat' array
                vert_index = lo*m + la + 1;
                text_index = lo*(m+1) + la + 1;

                rad_la = ((la/m)*2*pi)
                v = 1-(la / m);

                # Generate point on the torus at a given longitude/latitude
                gen_p = Vec3((-r * cos(rad_la) + 1) * sin(rad_lo), r * sin(rad_la), (-r * cos(rad_la) + 1) * cos(rad_lo));

                # Push vertices, texcoords, normals
                push!(positions, gen_p);
                push!(textures, Vec2(u, v));
                push!(normals, normalize(gen_p - inner_point));

                # Store 'flat' array indices for lookup
                vert_indices[lo+1, la+1] = vert_index;
                text_indices[lo+1, la+1] = text_index;
            end

            # Push texcoord that lies on the seam, 1 for each row
            push!(textures, Vec2(u, 0));
            text_indices[lo+1, m+1] = (lo+1)*(m+1);
        end
        
        # Push one last column of seam texcoords
        max_row = n*(m+1)+1;
        for g = 0:m
            push!(textures, Vec2(1,1-(g/m)));
            text_indices[n+1, g+1] = max_row + g;
        end
        
        # Generate all triangles
        for x = 1:n
            for y = 1:m
                # Wrapped so triangles can share vertices across the seam
                xp1 = (x+1 > n ? x + 1 - n : x + 1);
                yp1 = (y+1 > m ? y + 1 - m : y + 1);

                # Calculate indices for normals and vertices
                nv_indices = [vert_indices[xp1, y], vert_indices[x,y], vert_indices[xp1,yp1]];
                nv_indices_alt = [vert_indices[xp1, yp1], vert_indices[x,y], vert_indices[x,yp1]];

                # Calculate texcoord indices (different from the above as they do not wrap)
                t_indices = [text_indices[x+1, y], text_indices[x,y], text_indices[x+1,y+1]];
                t_indices_alt = [text_indices[x+1, y+1], text_indices[x,y], text_indices[x,y+1]];

                # Create and push triangles using the above indices
                push!(triangles, OBJTriangle(nv_indices, t_indices, nv_indices));
                push!(triangles, OBJTriangle(nv_indices_alt, t_indices_alt, nv_indices_alt));

            end
        end

        # Create and return the mesh
        return OBJMesh(positions, textures, normals, triangles);
    end

   #=
    = vectors_compare(indices1, vectors1, indices2, vectors2, epsilon)
    =
    = See whether two sets of indices point to the same set of vectors.
    =#
    function vectors_compare(indices1, vectors1, indices2, vectors2, epsilon)
        if(length(indices1) != length(indices2))
            return false;
        end
        
        # Check if every value value from the first index set is contained in the second index set.
        indices_count = length(indices1);
        
        for v1 = 1:indices_count
            # Track if the current vector was able to be matched
            vect_matched = false;

            for v2 = 1:indices_count
                # Compare if the two current vectors are equal (vary less than epsilon)
                if (norm(vectors1[indices1[v1]] - vectors2[indices2[v2]]) < epsilon)
                    vect_matched = true;
                    break;
                end
            end

            # If one vector could not be found, return false as the sets are not equal.
            if(!vect_matched)
                return false;
            end
        end

        return true;
    end

   #=
    = mesh_compare(mesh1, mesh2, verbose, epsilon)
    =
    = Compare whether two meshes are equal. Additionally, all errors can be logged and the precision of the comparison can be controlled.
    =#
    function mesh_compare(mesh1::OBJMesh, mesh2::OBJMesh, verbose::Bool = false, epsilon::Float64 = 0.001) :: Bool
        # Keep track of whether the mesh is valid.
        # This is important as while the mesh may be invalid early on, we want to print all violations.
        meshes_equal = true;

        # Check if the meshes have the same triangle count. If not, print an error message
        if(length(mesh1.triangles) != length(mesh2.triangles))
            if(verbose)
                verbose && println("[Compare] Mesh 1 and mesh 2 do not have the same amount of triangles.");
            end
            meshes_equal = false;
        end

        # Checks for whether each mesh has any uvs or normals
        has_uvs_1 = length(mesh1.uvs) > 0;
        has_uvs_2 = length(mesh2.uvs) > 0;
        has_norms_1 = length(mesh1.normals) > 0;
        has_norms_2 = length(mesh2.normals) > 0;

        # If only one mesh has UVs, error.
        if(has_uvs_1 != has_uvs_2)
            verbose && println("[Compare] Only one of the meshes has UVs.");
            meshes_equal = false
        end
        
        # If only one mesh has normals, error.
        if(has_norms_1 != has_norms_2)
            println("[Compare] Only one of the meshes has normals.");
            meshes_equal = false;
        end

        # Look for all triangles of mesh 1 in mesh 2, then repeat but reversed.
        for order = 1:2

            # Swap meshes for the second pass
            if(order == 2)
                buffer = mesh1;
                mesh1 = mesh2;
                mesh2 = buffer;
            end

            # Keep track of triangle number for logging purposes
            index = 0;
            
            # Check if each triangle in the current mesh is in the other mesh.
            for triangle1 in mesh1.triangles
                index += 1;

                # Track whether the current triangle has been matched.
                found_match = false;

                # Check current triangle against all other triangles in the other mesh.
                for triangle2 in mesh2.triangles
                    matched = true;

                    # Check the two current triangles' vertices against one another
                    if(length(triangle1.positions) == length(triangle2.positions))
                        matched =  vectors_compare(triangle1.positions, mesh1.positions, triangle2.positions, mesh2.positions, epsilon);
                    else
                        matched = false;
                    end
                    
                    if(!matched)
                        
                        continue;
                    end


                    # Check the two current triangles' UVs against one another
                    if(has_uvs_1 && has_uvs_2)
                        matched = vectors_compare(triangle1.uvs, mesh1.uvs, triangle2.uvs, mesh2.uvs, epsilon);
                    elseif (has_uvs_1 || has_uvs_2)
                        matched = false;
                    end

                    if(!matched)
                        continue;
                    end


                    # Check the two current triangles' normals against one another
                    if(has_norms_1 && has_norms_2)
                        matched = vectors_compare(triangle1.normals, mesh1.normals, triangle2.normals, mesh2.normals, epsilon);
                    elseif (has_norms_1 || has_norms_2)
                        matched = false;
                    end

                    if(!matched)
                        continue;
                    end

                    # If a match was found, indicate so and break.
                    if(matched)
                        found_match = true;
                        break;
                    end
                end

                # If a triangle was ever not found, lock valid to false.
                meshes_equal = meshes_equal && found_match;

                # Report if a triangle could not be found
                if(!found_match)
                    verbose && println("[Compare] Couldn't find a match for mesh " * string(order) * " triangle " * string(index) * " in mesh "*string(3 - order));
                end
            end
        end

        # Return whether the meshes are equal
        return meshes_equal;
    end

        """ 
    compress_obj(infile, outfile)
    Takes in an obj file, removing all duplicate vertices, texture coords, and normals.
    Updates indices as necessary.
    """
    function compress_obj(infile, outfile)
        mesh = read_obj(infile)

        unique_positions = Dict()
        new_position_indices = Dict()
        
        unique_normals = Dict()
        new_normal_indices = Dict()

        unique_uvs = Dict()
        new_uvs_indices = Dict()

        # iterate through the positions, if we find a new one, add it to the dictionary with its new index
        # if we find a position that we have seen before, add a new key to new_position_indices referring to the new index
        dup = 0
        for i=1:length(mesh.positions)
            p = mesh.positions[i]
            if !(p in keys(unique_positions))
                unique_positions[p] = length(unique_positions) + 1
            else
                new_position_indices[i] = unique_positions[p] + dup
                dup += 1
            end
        end


        # iterate through the normals, if we find a new one, add it to the dictionary with its new index
        # if we find a normal that we have seen before, add a key to new_normal_indices referring old index to the new index
        dup = 0
        for i=1:length(mesh.normals)
            n = mesh.normals[i]
            if !(n in keys(unique_normals))
                unique_normals[n] = length(unique_normals) + 1
            else
                new_normal_indices[i] = unique_normals[n] + dup
                dup += 1
            end
        end

        # iterate through the uvs, if we find a new one, add it to the dictionary with its new index
        # if we find a uvs that we have seen before, add a new key to new_uvs_indices referring to the new index
        dup = 0
        for i=1:length(mesh.uvs)
            uv = mesh.uvs[i]
            if !(uv in keys(unique_uvs))
                unique_uvs[uv] = length(unique_uvs) + 1
            else
                new_uvs_indices[i] = unique_uvs[uv] + dup
                dup += 1
            end
        end


        for i = 1:length(mesh.triangles)
            # check if vertex indices are outdated, update if needed
            triangle = mesh.triangles[i]
            for j = 1:length(triangle.positions)
                old_index = triangle.positions[j]
                if old_index in keys(new_position_indices)
                    mesh.triangles[i].positions[j] = new_position_indices[old_index]
                end
                #update indices that are affected by removing other points
                count = 0
                for removed_index in keys(new_position_indices)
                    if mesh.triangles[i].positions[j] > removed_index
                        count += 1
                    end
                end
                mesh.triangles[i].positions[j] -= count
            end

            #check if normal indices are outdated, update if needed
            for j = 1:length(triangle.normals)
                old_index = triangle.normals[j]
                if old_index in keys(new_normal_indices)
                    mesh.triangles[i].normals[j] = new_normal_indices[old_index]
                end
                #update indices that are affected by removing other points
                count = 0
                for removed_index in keys(new_normal_indices)
                    if mesh.triangles[i].normals[j] > removed_index
                        count += 1                    
                    end
                end
                mesh.triangles[i].normals[j] -= count
            end

            #check if uv indices are outdated, update if needed
            for j = 1:length(triangle.uvs)
                old_index = triangle.uvs[j]
                
                if old_index in keys(new_uvs_indices)
                    mesh.triangles[i].uvs[j] = new_uvs_indices[old_index]
                end
                #update indices that are affected by removing other points
                count = 0
                for removed_index in keys(new_uvs_indices)
                    if mesh.triangles[i].uvs[j] > removed_index
                        count += 1                    
                    end
                end
                mesh.triangles[i].uvs[j] -= count
            end

        end

        # delete old duplicate coordinates
        new_positions = []
        new_normals = []
        new_uvs = []

        for i = 1:length(mesh.positions)
            if !(i in keys(new_position_indices))
                push!(new_positions, mesh.positions[i])
            end
        end

        for i = 1:length(mesh.normals)
            if !(i in keys(new_normal_indices))
                push!(new_normals, mesh.normals[i])
            end
        end

        for i = 1:length(mesh.uvs)
            if !(i in keys(new_uvs_indices))
                push!(new_uvs, mesh.uvs[i])
            end
        end

        mesh.uvs = new_uvs
        mesh.positions = new_positions
        mesh.normals = new_normals

        write_obj(outfile, mesh)
    end



   #= estimate_normals(mesh::OBJMesh)
    = 
    = Estimates normals for a given mesh. Overwrites any existing normals and returns a new OBJMesh object.
    =#
    function estimate_normals(mesh::OBJMesh)
        # Create mesh data arrays
        normals = zeros(Vec3, size(mesh.positions)[1]);
        positions = mesh.positions;
        textures = mesh.uvs;
        triangles = [];

        for triangle in mesh.triangles
            # Clear triangle normals, use vertex indices instead as per-vertex normals are the target
            triangle.normals = triangle.positions;
            indices = triangle.normals;

            # Calculate the triangle's normal usng the cross product
            calc_normal = normalize(cross(positions[indices[1]] - positions[indices[3]], positions[indices[2]] - positions[indices[3]]));

            # Add the normal to all vertices in the triangle
            for n in indices
                normals[n] += calc_normal;
            end

            # Add the modified triangle to the new array
            push!(triangles, triangle);
        end

        # Normalize all normal vectors as they are probably huge
        len_normals = size(normals)[1];
        for i = 1:len_normals
            normals[i] = normalize(normals[i]);
        end

        # Create and return a new mesh
        return OBJMesh(positions, textures, normals, triangles);
    end

end # module WWUMeshes


