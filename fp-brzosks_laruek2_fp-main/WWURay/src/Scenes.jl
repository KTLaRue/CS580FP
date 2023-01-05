module Scenes

export HitRecord, Sphere, Scene, TriangleMesh, ray_intersect, create_triangles, Cloud
#export has_uvs, has_normals, get_vertex, get_uv, get_normal

using LinearAlgebra
#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..WWUMeshes
using ..Materials



#####################################
###### Generic Scene Data Type ######
#####################################
struct Scene
    background::RGB{Float32}
    objects::Array{Any,1}
    lights::Array{Any,1}
    air::Array{Any,1}
end

""" Structure to store data about an intersection
of a ray with an object (a "hit")."""
mutable struct HitRecord
    t::Float64
    intersection::Vec3
    normal::Vec3
    uv::Union{Vec2,Nothing}
    object
end

# Abstract ray-object intersection function:
# Needs to be implemented for each type of object to be rendered
""" Intersect a ray with an object.
Returns a HitRecord with info about the intersetion, or nothing if
the ray doesn't intersect with the object. """
function ray_intersect(ray::Ray, object) end


##################
##### Sphere #####
##################

# Data type:
struct Sphere
    center::Vec3
    radius::Float64
    material::Material
end

struct Cloud
    center::Vec3
    radius::Float64
    density::Float64
end

""" Ray-sphere intersection. """
function ray_intersect(ray::Ray, object::Sphere)

    # find the A, B, C for solving the quadratic formula
    A = dot(ray.direction, ray.direction)
    B = 2 * (dot(ray.origin - object.center, ray.direction))
    C = dot(ray.origin - object.center, ray.origin - object.center) - object.radius ^ 2

    # calculate the discriminant
    discrim = B ^ 2 - 4 * A * C

    # if < 0 we can only have imaginary solutions -- ray does not intersect sphere
    if discrim < 0
        return nothing
    # if == 0 then there can only be one solution -- ray intersects sphere at ony point
    elseif discrim == 0
        t = -1 * B / (2 * A)
        if t < 0
            return nothing
        end
    # if > 0 then there are 2 solutions -- ray intersects at 2 points
    else
        # calculate the 2 t possible values
        t_plus = (-1 * B + sqrt(discrim)) / (2 * A)
        t_minus = (-1 * B - sqrt(discrim)) / (2 * A)
        # if both are negative, return nothing -- behind the ray's origin
        if t_plus < 0 && t_minus < 0
            return nothing
        else
            # calculate the min t that is greater than 0
            t_min = min(t_plus, t_minus)
            if t_min < 0
                t = max(t_plus, t_minus)
            else
                t = t_min
            end
        end   
    end


    # calculate the intersection and normal for the ray-impact point  
    intersection = ray.origin + (t * ray.direction)
    normal = normalize(intersection - object.center)

    #adjust for sphere not centered at the origin
    x = intersection[1] - object.center[1]
    y = intersection[2] - object.center[2]
    z = intersection[3] - object.center[3]
    
    # calculate longitude and latitude
    lon = -atan(z,x) + pi/2
    lat = asin(y / object.radius)

    #calculate u, v from lat/lon 
    u = lon / (2 * pi) + 0.5
    v = (lat / pi) + 0.5
    uvs = Vec2(u, v)

    hit_info = HitRecord(t, intersection, normal, uvs, object)

    return hit_info

    ################################################################
    # TODO 9c - modify above to fill in Hitrec's texture coordinates
    ################################################################
end


###########################
###### Triangle Mesh ######
###########################

""" Data type: stores the OBJTriangle, a reference to its Mesh
object, and the material it should be rendered with. """
struct Triangle
    geometry::OBJTriangle
    mesh::OBJMesh
    material
end

""" Return an Array of all Triangles belonging to the given mesh, assigning
each one the given material. """
function create_triangles(mesh::OBJMesh, material)
    [Triangle(f, mesh, material) for f in mesh.triangles]
end

""" Some helper functions that make for easier access to triangle data: """
function get_vertex(tri::Triangle, i)
    tri.mesh.positions[tri.geometry.positions[i]]
end

function has_uvs(tri::Triangle)
    length(tri.geometry.uvs) == 3
end

function get_uv(tri::Triangle, i)
    tri.mesh.uvs[tri.geometry.uvs[i]]
end

function has_normals(tri::Triangle)
    length(tri.geometry.normals) == 3
end

function get_normal(tri::Triangle, i)
    tri.mesh.normals[tri.geometry.normals[i]]
end

function ray_intersect(ray::Ray, object::Triangle)
    
    ##########
    # TODO 7 #
    ##########
    # Your implementation:
    #
    # solve for beta, gamma, and t using Cramer's Rule:
    pointA = object.mesh.positions[object.geometry.positions[1]]
    pointB = object.mesh.positions[object.geometry.positions[2]] 
    pointC = object.mesh.positions[object.geometry.positions[3]] 
    dir = ray.direction
    eye = ray.origin
    
    a = pointA[1] - pointB[1]
    b = pointA[2] - pointB[2]
    c = pointA[3] - pointB[3]
    d = pointA[1] - pointC[1]
    e = pointA[2] - pointC[2]
    f = pointA[3] - pointC[3]
    g = dir[1]
    h = dir[2]
    i = dir[3]
    j = pointA[1] - eye[1]
    k = pointA[2] - eye[2]
    l = pointA[3] - eye[3]
    M = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g)

    beta = (j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g)) / M
    gamma = (i * (a * k - j * b) + h * (j * c - a * l) + g * (b * l - k * c)) / M
    alpha = 1 - beta - gamma
    t = -(f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c)) / M

    # make sure that we have hit a point that is within the triangle interior or on the edge.
    if beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && alpha <= 1 && alpha >= 0
        intersection = eye + (t * dir)
        uv = nothing
        # calculate the interpolated uv from the barycentric coords
        if length(object.geometry.uvs) != 0
            uvA = object.mesh.uvs[object.geometry.uvs[1]]
            uvB = object.mesh.uvs[object.geometry.uvs[2]]
            uvC = object.mesh.uvs[object.geometry.uvs[3]]
            uv = alpha * uvA + beta * uvB + gamma * uvC
        end
        if length(object.geometry.normals) == 0
            # calculate the normal for a mesh with no normal values stored at vertices
            side_one = pointB - pointA
            side_two = pointC - pointA
            normal = normalize(cross(side_one, side_two))
        else
            # calculate interpolated normal from barycentric coordinates
            normA = object.mesh.normals[object.geometry.normals[1]]
            normB = object.mesh.normals[object.geometry.normals[2]]
            normC = object.mesh.normals[object.geometry.normals[3]]
            normal = normalize(alpha * normA + beta * normB + gamma * normC)
        end

        hit_info = HitRecord(t, intersection, normal, uv, object)

        return hit_info 
    else
        return nothing
    end

    ##############
    # END TODO 7 #
    ##############

    ################################################################
    # TODO 9c - modify above to fill in Hitrec's texture coordinates
    ################################################################
end

end # module Scenes
