""" 
Main module for CS480/580 A2 raytracer. Contains core raytracing algrithm,
while referencing several other modules that encapsulate supporting
functionality.
"""

module WWURay

export main

using FileIO
using Images
using StaticArrays
using LinearAlgebra

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
include("Lights.jl")
include("Materials.jl")
include("WWUMeshes.jl")
include("Scenes.jl")
include("Cameras.jl")
include("TestScenes.jl")

using .GfxBase
using .Lights
using .Materials

import .Scenes
import .Scenes.Scene
import .Scenes.HitRecord
import .Cameras
import .TestScenes

# Ray-Scene intersection:
""" Find the closest intersection point among all objects in the scene
along a ray, constraining the search to values of t between tmin and tmax. """
function closest_intersect(objects::Array{Any, 1}, ray::Ray, tmin, tmax)
    # call ray_intersect on each object in the scene
    # record HitRecord with the lowest t where t_min <= t <= t_max
    # return the HitRecord with corresponding t
    closest_t = Inf
    closest_hit = nothing
    # iterate over all objects in scene
    for obj in objects
        # find the (possible) intersection between ray and obj
        hit = Scenes.ray_intersect(ray, obj)
        # check if we have found a new min dist with t in appropriate range -- update accordingly
        if hit != nothing && hit.t > tmin && hit.t < tmax && hit.t < closest_t
            closest_t = hit.t
            closest_hit = hit
        end
    end

    # if no objects were intersected -- closest_hit == nothing
    return closest_hit
end

""" find_point_cyl(o, t, l)
Given an 3d origin o, angle t in rads, and length l, find the coordinates of the end of the
vector originating from o
"""
function find_point_hemi(t, l)
    x = cos(t)*sin(l)
    y = cos(l)
    z = sin(t)*sin(l)
    return Vec3(x, y, z)
end

"""samples random directions based on unit circle around intersection for diffuse color"""
function diffuse_sample(scene::Scene, hitrec::HitRecord, count::Int, view::Vec3, depth::Int, maxDepth::Int)
    #How many bounces wil we have max
    if depth > maxDepth
        return RGB{Float32}(0, 0, 0)
    end

    point = hitrec.intersection
    normal = hitrec.normal
    imaterial = hitrec.object.material
    ishader = imaterial.shading_model

    ray_color = RGB{Float32}(0.0, 0.0, 0.0)

    #Change Basis
    modified_view = nothing
    if view[1] >= view[2]
        modified_view = Vec3(0.0, view[2], view[3])
    elseif view[2] >= view[3]
        modified_view = Vec3(view[1], 0.0, view[3])
    else
        modified_view = Vec3(view[1], view[2], 0.0)
    end  

    v = normalize(normal)
    u = normalize(cross(v, normalize(modified_view)))
    w = cross(u,v)

    for i in 1:count
        #get random point in unit circle around intersection on object plane
        xz = rand(Float64,1)[1] * pi * 2
        y = rand(Float64,1)[1] * pi/2
        direct_point = find_point_hemi(xz, y)

        #project point to unit sphere
        generated_vec = direct_point[1] * u + direct_point[2] * v + direct_point[3] * w

        #geberate ray and find if it interesects anything
        generated_ray = Ray(point, generated_vec)
        closest_hitrec = closest_intersect(scene.objects, generated_ray, 1e-8, 40)
        if closest_hitrec != nothing
            for light in scene.lights
                int_point = closest_hitrec.intersection
                int_norm = closest_hitrec.normal
                material = closest_hitrec.object.material
                shader = material.shading_model

                #Consider the surface as a light source
                local_color = determine_color(shader, material, generated_ray, closest_hitrec, scene)
                if local_color != RGB{Float32}(0.0, 0.0, 0.0)
                    lit = DirectionalLight(local_color, generated_ray.direction)

                    ray_color += shade_light(ishader, imaterial, generated_ray, hitrec, lit, scene)
                    ray_color += diffuse_sample(scene, closest_hitrec, count, generated_ray.direction, depth+1, maxDepth)
                end
            end
        end
    end

    #scale color by how many rays we generate
    ray_color = ray_color/count
                    
    return ray_color
end

"""colects random direction in unit circle around intersection with more probability around mirror direction for glossy surfaces"""
function glossy_sample(scene::Scene, hitrec::HitRecord, count::Int, view::Vec3, gloss::Float64, depth::Int, maxDepth::Int)
    #How many bounces wil we have max
    if depth > maxDepth
        return RGB{Float32}(0, 0, 0)
    end

    point = hitrec.intersection
    normal = hitrec.normal
    imaterial = hitrec.object.material
    ishader = imaterial.shading_model

    ray_color = RGB{Float32}(0.0, 0.0, 0.0)

    #Determine how many perferred rays vs. nonperferred rays to generate
    #The higher the mirror coefficient, the more rays in the reflected direction
    pCount = round(Int, count*gloss)
    nCount = count - pCount

    #Find mirrored Ray
    D = view * -1
    reflect = normalize(2 * normal * dot(normal, D) - D)

    modified_norm = nothing
    if normal[1] >= normal[2]
        modified_norm = Vec3(0.0, normal[2], normal[3])
    elseif normal[2] >= normal[3]
        modified_norm = Vec3(normal[1], 0.0, normal[3])
    else
        modified_norm = Vec3(normal[1], normal[2], 0.0)
    end  

    #Basis Change
    v = normalize(reflect)
    u = normalize(cross(v, normalize(normal)))
    w = cross(u,v)

    #Create symetric perfered rays list
    new_ray = Ray(point, reflect)
    ang = 1-gloss

    perfered_set = []
    for i in 1:pCount
        #get random point in unit circle around intersection on object plane
        xz = rand(Float64,1)[1] * pi * 2
        y = (rand(Float64,1)[1]*ang) * pi/2
        direct_point = find_point_hemi(xz, y)

        generated_vec = direct_point[1] * u + direct_point[2] * v + direct_point[3] * w

        #generate ray and find if it interesects anything
        generated_ray = Ray(point, generated_vec)
        if dot(normal, generated_vec) > 0.01
            push!(perfered_set, generated_ray)
        else
            i = i-1
        end
    end

    v = normalize(normal)
    u = normalize(cross(v, normalize(view)))
    w = cross(u,v)

    nonperfered_set = []
    for i in 1:nCount
        #get random point in unit circle around intersection on object plane
        xz = rand(Float64,1)[1] * pi * 2
        y = rand(Float64,1)[1] * pi/2
        direct_point = find_point_hemi(xz, y)

        #project point to unit sphere
        generated_vec = direct_point[1] * u + direct_point[2] * v + direct_point[3] * w
        #generate ray and find if it interesects anything
        generated_ray = Ray(point, generated_vec)
        if dot(normal, generated_vec) > 0.01
            push!(nonperfered_set, generated_ray)
        else
            i = i-1
        end
    end

    if pCount > 0
        for r in perfered_set
            closest_hitrec = closest_intersect(scene.objects, r, 1e-8, 40)
            if closest_hitrec != nothing
                for light in scene.lights
                    int_point = closest_hitrec.intersection
                    int_norm = closest_hitrec.normal
                    material = closest_hitrec.object.material
                    shader = material.shading_model

                    #Consider the surface as a light source
                    local_color = determine_color(shader, material, r, closest_hitrec, scene)
                    if local_color != RGB{Float32}(0.0, 0.0, 0.0)
                        lit = DirectionalLight(local_color, r.direction)

                        ray_color += shade_light(ishader, imaterial, r, hitrec, lit, scene)*(nCount/pCount)
                        ray_color += glossy_sample(scene, closest_hitrec, count, r.direction, gloss, depth+1, maxDepth)
                    end
                end
            end
        end
    end
    if nCount > 0
        for r in nonperfered_set
            closest_hitrec = closest_intersect(scene.objects, r, 1e-8, 40)
            if closest_hitrec != nothing
                for light in scene.lights
                    int_point = closest_hitrec.intersection
                    int_norm = closest_hitrec.normal
                    object = closest_hitrec.object
                    material = object.material
                    shader = material.shading_model

                    #Consider the surface as a light source
                    local_color = determine_color(shader, material, r, closest_hitrec, scene)
                    if local_color != RGB{Float32}(0.0, 0.0, 0.0)
                        lit = DirectionalLight(local_color, r.direction)

                        ray_color += shade_light(ishader, imaterial, r, hitrec, lit, scene)*(pCount/nCount)
                        ray_color += glossy_sample(scene, closest_hitrec, count, r.direction, gloss, depth+1, maxDepth)
                    end
                end
            end
        end
    end
    
    ray_color = ray_color/count           
    return ray_color
end

""" Trace a ray from orig along ray through scene, using Whitted recursive raytracing 
limited to rec_depth recursive calls. """
function traceray(scene::Scene, ray::Ray, tmin, tmax, volumetric::Bool, global_fog::Bool, vCount::Int, 
                    illum::Bool, gCount::Int, gDepth::Int, rec_depth = 1)

    closest_hitrec = closest_intersect(scene.objects, ray, tmin, tmax)

    # if we dont hit anything, dont return a color
    if closest_hitrec == nothing
        return scene.background
    end

    object = closest_hitrec.object
    P = closest_hitrec.intersection
    normal = closest_hitrec.normal
    material = object.material
    shader = material.shading_model
    light_ray = P - ray.origin

    local_color = determine_color(shader, material, ray, closest_hitrec, scene)

    #if volumetric is true
    if volumetric
        if vCount == 0
            vCount = 1
        end
        steps = vCount

        #add volumetric lighting
        light_seg_step = light_ray / steps

        mag_seg = length(light_seg_step)

        #need to check if this should be somthing else
        extinctionCoef = .1
        ambientLight = 1
        scatteringCoef = .1 #refelects almost density of what is in the air of the room
        inScattering = 0.0
        outScattering = 0.0
        currentLight = 0.0
        illumination = RGB{Float32}(0.0, 0.0, 0.0)
        transmittance = 1

        sample_pos = ray.origin
        for segment in 1:steps
            #get our current location
            sample_pos += light_seg_step 

            for light in scene.lights
                is_direct = direct_light(light, scene, sample_pos)
                # cur_density = find_density(scene, sample_pos, global_fog)
                # transmittance = transmittance * exp(-cur_density * extinctionCoef * mag_seg)

                #if our point is directly in a light path
                if is_direct == nothing
                    #get values for light scattering
                    cur_density = find_density(scene, sample_pos, global_fog)
                    transmittance = transmittance * exp(-cur_density * extinctionCoef * mag_seg)

                    #Needs check for directional vs point light
                    inScattering = light.intensity
                    if isa(light, WWURay.Lights.PointLight)
                        inScattering = (light.intensity)/(0.2*norm(light.position-sample_pos))
                    end

                    outScattering = scatteringCoef * cur_density
                    currentLight = inScattering * outScattering
                    illumination += transmittance * currentLight * mag_seg
                end
            end
        end
        local_color = local_color + illumination
    end

    #if we have global illumination
    if illum
        if gCount == 0
            gCount = 1
        end
        glossy_coeff = closest_hitrec.object.material.glossy_coeff
        if glossy_coeff == 0.0
            #if we have a diffuse surface
            global_color = diffuse_sample(scene, closest_hitrec, gCount, light_ray, 0, gDepth)
            local_color = local_color + global_color
        elseif (glossy_coeff > 0.0)
            #if we have glossy surface
            global_color = glossy_sample(scene, closest_hitrec, gCount, light_ray, glossy_coeff, 0, gDepth)
            local_color = local_color + global_color
        end
    end

    # return color if we have hit max reflection depth or if there is no mirroring on the object
    if rec_depth == 8 || material.mirror_coeff == 0
        return local_color 
    end

    #if we have a mirror and/or bouncing rays
    D = ray.direction * -1
    N = normal
    reflect = 2 * N * dot(N, D) - D
    new_ray = Ray(P, reflect)
    
    #recursively calculate color from new reflection ray
    reflected_color = traceray(scene, new_ray, 1e-8, tmax, volumetric, global_fog, vCount, illum, gCount, rec_depth + 1)

    return material.mirror_coeff * reflected_color + (1 - material.mirror_coeff) * local_color
end

"""finds the density of the room at a given 3D point based on sections of localized fog/mist/smoke"""
function find_density(scene::Scene, pos::Vec3, global_fog::Bool)
    #generally, air has a density so small it could be 0
    density = 0
    for section in scene.air
        dist = norm(pos - section.center)
        if dist <= section.radius
            density += section.density
        end
    end

    #if you want a global density of .1
    if global_fog
        density = density + .1
    end
        
    return density
end

"""returns nothing if point is in the direct path of a given light, else it returns the hitrec of the object blocking"""
function direct_light(lights::Light, scene::Scene, point::Vec3)
    result = nothing
    if isa(lights, WWURay.Lights.PointLight)
        ray = Ray(point, lights.position - point)
        result = closest_intersect(scene.objects, ray, 1e-8, 1)
    else
        ray = Ray(point, lights.direction)
        result = closest_intersect(scene.objects, ray, 1e-8, Inf)
    end
    return result
end

""" Determine the color of interesction point described by hitrec 
Flat shading - just color the pixel the material's diffuse color """
function determine_color(shader::Flat, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    get_diffuse(material, hitrec.uv)
end
""" Normal shading - color-code pixels according to their normals """
function determine_color(shader::Normal, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    normal_color = normalize(hitrec.normal) / 2 .+ 0.5
    RGB{Float32}(normal_color...)
end

""" Determine the color of a physical (Lambertian, BlinnPhong, etc.) surface """
function determine_color(shader::PhysicalShadingModel, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)

    # start with a black color value
    # for each light in the scene:
    #   determine the light's contribution (by calling shade_light)
    #   add the light's contribution into the color
    # return the resulting color

    color = RGB{Float32}(0.0, 0.0, 0.0)
    for light in scene.lights
        if (!is_shadowed(scene, light, hitrec.intersection))
            color += shade_light(shader, material, ray, hitrec, light, scene)
        end
    end
    return color
end

""" shade_light(shader, material, ray, hitrec, light, scene)
Determine the color contribution of the given light along the given ray.
Color depends on the material, the shading model (shader), properties of the intersection 
given in hitrec, """
function shade_light(shader::Lambertian, material::Material, ray::Ray, hitrec, light, scene)

    #calculate color contribution using the lambertian formula:
    k_d = get_diffuse(material, hitrec.uv)
    I = light.intensity 
    color = RGB{Float32}(k_d.r * I.r, k_d.g * I.g, k_d.b * I.b)
    L = color * max(0, dot(hitrec.normal, normalize(light_direction(light, hitrec.intersection))))
    return L
end

""" Blinn-Phong surface shading """
function shade_light(shader::BlinnPhong, material::Material, ray::Ray, hitrec, light, scene)

    k_d = get_diffuse(material, hitrec.uv)
    k_s = shader.specular_color
    n = hitrec.normal
    I = light.intensity
    color = RGB{Float32}(k_d.r * I.r, k_d.g * I.g, k_d.b * I.b) 
    l = normalize(light_direction(light, hitrec.intersection))
    v = -normalize(ray.direction)
    h = normalize(v + l)
    p = shader.specular_exp

    #calculate color contribution using the blinn phong formula:
    L = color * max(0, dot(n, l)) + k_s * I * max(0, dot(n, h))^p

    return L
end


""" Determine whether point is in shadow wrt light """
function is_shadowed(scene, light::DirectionalLight, point::Vec3)

    ray = Ray(point, normalize(light.direction))
    hit = closest_intersect(scene.objects, ray, 1e-8, Inf)
    # if we hit an object, then we know there is something obscuring the light source
    if hit != nothing
        return true
    else    
        return false
    end
end

function is_shadowed(scene, light::PointLight, point::Vec3)
    ray = Ray(point, light.position - point) 
    hit = closest_intersect(scene.objects, ray, 1e-8, 1)
    # if we hit an object, then we know there is something obscuring the light source
    if hit != nothing
        return true
    else    
        return false
    end
end

# Main loop:
function main(scene, camera, height, width, outfile, volumetric, global_fog, local_fog, vCount, illum, gCount, gDepth)
    """parameter notes
    the folowing are boolean values that denote curtain features
    volumetric  - is volumetric lighting on (T) or off (F)
    global_fog  - do you want the scene to have a global fog? yes (T) no (F)
    local_fog   - is there pre determined localized fog? yes (T) no (F)
    vCount      - sampling rate for volumetric lighting
    illum       - do you want global illumination?
    gCount      - sampling rate for global illumination
    gDepth      - number of bounces in global illumination
    """

    # get the requested scene and camera
    scene = TestScenes.get_scene(scene, local_fog)
    camera = TestScenes.get_camera(camera, height, width)

    # Create a blank canvas to store the image:
    canvas = zeros(RGB{Float32}, height, width)

    # Pseudocode:
    #   loop over all pixels in the image
    #   for each pixel, get a viewing ray from the camera
    #   then call traceray to determine its color
    #
    h, w = size(canvas)
    for i=1:h
        for j=1:w
            viewing_ray = Cameras.pixel_to_ray(camera, i, j)
            color = traceray(scene, viewing_ray, 1, Inf, volumetric, global_fog, vCount, illum, gCount, gDepth)
            canvas[i, j] = color
        end
    end

    # clamp canvas to valid range:
    clamp01!(canvas)
    save(outfile, canvas)
end

end # module WWURay

