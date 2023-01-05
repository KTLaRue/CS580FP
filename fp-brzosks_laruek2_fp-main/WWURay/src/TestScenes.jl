module TestScenes
#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..Scenes
using ..Materials
using ..Lights
using ..WWUMeshes
using ..Cameras

# helpful things:
make_diffuse(color) = Material(Lambertian(), 0.0, nothing, color)
black = RGB{Float32}(0,0,0)
red = RGB{Float32}(1,0,0)
green = RGB{Float32}(0,1,0)
blue = RGB{Float32}(0,0,1)
white = RGB{Float32}(1,1,1)
purple = RGB{Float32}(1,0,1)

function camera_1(img_height, img_width)
    CanonicalCamera(img_height, img_width)
end

function camera_2(img_height, img_width)
    eye = Vec3(20, 4, 10)
    view = Vec3(-1, 0, -5) - eye
    up = Vec3(0, 1, 0)
    focal = 8.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end

function camera_3(img_height, img_width)

    Cameras.PerspectiveCamera(
                 Vec3(-1, 0.8, -1.2),  # eye::Vec3
                 Vec3(1, -1, -1), # view::Vec3
                 Vec3(0, 1, 0),   # up::Vec3
                 0.3,     # focal::Real
                 img_height, # canv_height::Int
                 img_width) # canv_width::Int)
end

cameras = [camera_1, camera_2, camera_3]

function get_camera(i, img_height, img_width)
    cameras[i](img_height, img_width)
end

function get_scene(i)
    scene[i]()
end

function get_scene(i, local_fog)
    scenes[i](local_fog)

end

function scene_1()
    bg = RGB{Float32}(0.95, 0.95, 0.95)
    objs = [Sphere(Vec3(0, 0, -5), 1, Material(Flat(), 0.0, nothing, RGB{Float32}(0.73,0,0.17)))]
    lights = [PointLight(0.8, Vec3(0,0,0))]
    air=[]
    Scene(bg, objs, lights, air)
end

function scene_2()
    bg = black
    objs = [
            Sphere(Vec3( 2, 0, -8), 1, Material(Lambertian(), 0.0, nothing, white)),
            Sphere(Vec3(-2, 0, -8), 2, Material(Lambertian(), 0.0, nothing, blue))
           ]

    lights = [ DirectionalLight(1.0, Vec3(1, 0.5, -0.1)) ]
    air=[]
    Scene(bg, objs, lights, air)
end

function scene_3()
    bg = black
    mat = Material(Lambertian(), 0.0, nothing, white)
    objs = [
            Sphere(Vec3( -2, 1, -8), 1, mat),
            Sphere(Vec3(  2, 1, -8), 1, mat)
           ]

    lights = [ PointLight(1.0, Vec3(0, 5, -8.5)) ]

    air=[]
    Scene(bg, objs, lights, air)
end

function scene_4()
    bg = black
    mat1 = Material(BlinnPhong(white, 10), 0.0, nothing, white)
    mat2 = Material(BlinnPhong(white, 10), 0.0, nothing, blue)
    mat3 = Material(BlinnPhong(red, 100), 0.0, nothing, blue)
    objs = [
            Sphere(Vec3( -2, -1, -8), 1, mat1),
            Sphere(Vec3( -1, 1, -8), 1, mat2),
            Sphere(Vec3(  0, -1, -8), 1, mat3),
            Sphere(Vec3(  1, 1, -8), 1, mat2),
            Sphere(Vec3(  2, -1, -8), 1, mat1),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.0, nothing, white))
           ]

    lights = [ PointLight(0.8, Vec3(0, 4, -8)),
               PointLight(0.2, Vec3(0, 0, 0)) ]

    air=[]
    Scene(bg, objs, lights, air)end

function scene_5()
    bg = black

    mat = Material(Lambertian(), 0.0, nothing, white)

    objs = [
            Sphere(Vec3( -1, 0, -6), 0.5, mat),
            Sphere(Vec3(  1, 0, -5), 0.5, Material(Lambertian(), 0.4, nothing, white)),
            Sphere(Vec3( -1, 0, -4), 0.5, mat),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, white)) # ground
           ]

    lights = [ DirectionalLight(0.6, Vec3(1, 1, 0)),
               PointLight(0.4, Vec3(0, 0, 0)) ]

    air=[]
    Scene(bg, objs, lights, air)end

function scene_6()
    bg = black

    r = Material(BlinnPhong(white, 10), 0.0, nothing, red)
    g = Material(BlinnPhong(white, 10), 0.0, nothing, green)
    b = Material(BlinnPhong(white, 10), 0.0, nothing, blue)
    refl = Material(Lambertian(), 0.6, nothing, white)

    objs = [
            #Sphere(Vec3(-10, 0, -1), 9.2, refl),
            Sphere(Vec3(-1,  -1.1, -3), 0.5, r),
            Sphere(Vec3( -0.5,  -1.0, -4), 0.5, g),
            Sphere(Vec3( 0,  -0.9, -5), 0.5, b),
            Sphere(Vec3( 5,  -1, -4), 4, refl),
            #Sphere(Vec3( 10,  0.1 , -1), 9.2, refl),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, white)) # floor
           ]

    lights = [ PointLight(0.6, Vec3(1, 10, -4)),
               PointLight(0.4, Vec3(0, 0, 0)) ]

    air=[]
    Scene(bg, objs, lights, air)end

""" Take the OBJMesh mesh and return an array of Triangles from the mesh
with the given material, after scaling the mesh positions by scale and moving
them by translation """
function mesh_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        mesh.positions[i] = mesh.positions[i] * scale + translation
    end

    create_triangles(mesh, material)
end

function scene_7()
    bg = black
    objs = []

    # add a bunny:
    bunny_mat = Material(Lambertian(), 0.0, nothing, RGB{Float32}(0.6, 0.5, 0.5))
    bunny = read_obj("data/bunny.obj")
    append!(objs, mesh_helper(bunny, bunny_mat, 1.0, Vec3(0.2, 0, -5)))

    # add a cube
    cube_mat = Material(Lambertian(), 0.6, nothing, white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 10.0, Vec3(-11.2, 0, 0)))

    lights = [ PointLight(0.5, Vec3(1,2,-5)),
               DirectionalLight(0.3, Vec3(0,0,1)),
               DirectionalLight(0.3, Vec3(0,1,1)),
               DirectionalLight(0.3, Vec3(1,1,1)),
               DirectionalLight(0.3, Vec3(0,1,0)) ]

    air=[]
    Scene(bg, objs, lights, air)
end

scene_8 = scene_7

function scene_9()
    bg = black

    objs = []

    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.2, nothing, RGB{Float32}(0.8, 0.8, 1.0))))

    sphere_material = Material(Lambertian(), 0.0, Texture("data/earth.png", false), nothing)
    push!(objs, Sphere(Vec3(-1.25, 0, -6), 1, sphere_material))

    sphere_m = sphere_mesh(32, 16)
    scale = 1.0
    translation = Vec3(1.25, 0, -6)
    for i in 1:length(sphere_m.positions)
        sphere_m.positions[i] = sphere_m.positions[i] * scale + translation
    end
    append!(objs, create_triangles(sphere_m, sphere_material))

    cube_mat = Material(Lambertian(), 0.0, Texture("data/1.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.5, Vec3(-1, -1, -3)))

    lights = [ DirectionalLight(0.4, Vec3(0,1,0)),
               DirectionalLight(0.8, Vec3(0.4,0.4,1)) ]

    air=[]
    Scene(bg, objs, lights, air)
end

function scene_10()
    bg = black
    objs = []

    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.4, nothing, RGB{Float32}(0.2, 0.8, .5))))

    mat1 = Material(Lambertian(), 0.0, Texture("data/opt7.png", false), nothing)
    push!(objs, Sphere(Vec3(-2.45, 0, -6), .3, mat1))

    mat2 = Material(Lambertian(), 0.0, Texture("data/opt1.png", false), nothing)
    push!(objs, Sphere(Vec3(-2.1, 1.3, -6), .3, mat2))

    mat3 = Material(Lambertian(), 0.0, Texture("data/opt2.png", false), nothing)
    push!(objs, Sphere(Vec3(-1.3, 2.1, -6), .3, mat3))
    
    
    mat4 = Material(Lambertian(), 0.0, Texture("data/opt3.png", false), nothing)
    push!(objs, Sphere(Vec3(0, 2.45, -6), .3, mat4))


    mat5 = Material(Lambertian(), 0.0, Texture("data/opt4.png", false), nothing)
    push!(objs, Sphere(Vec3(1.3, 2.1, -6), .3, mat5))

    mat6 = Material(Lambertian(), 0.0, Texture("data/opt5.png", false), nothing)
    push!(objs, Sphere(Vec3(2.1, 1.3, -6), .3, mat6))


    mat7 = Material(Lambertian(), 0.0, Texture("data/opt8.png", false), nothing)
    push!(objs, Sphere(Vec3(2.45, 0, -6), .3, mat7))

    #for placement
    sphere_material = Material(Lambertian(), 0.0, Texture("data/img1.png", false), nothing)
    push!(objs, Sphere(Vec3(0, 0, -15), 5, sphere_material))


    cube_mat = Material(Lambertian(), 0.0, Texture("data/opt6.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.4, Vec3(-.5, -1, -3)))

    lights = [ DirectionalLight(0.4, Vec3(0,1,0)),
               DirectionalLight(0.8, Vec3(0.4,0.4,1)) ]

    air=[]
    Scene(bg, objs, lights, air)

end

function scene_11()
    bg = black

    objs = []
 
    sun_material = Material(Lambertian(), 0.0, Texture("data/sun.jpg", false), nothing)
    push!(objs, Sphere(Vec3(-5, -1, -20), 10, sun_material))

    offset = Vec3(0, 0, 0)

    mercury_material = Material(Lambertian(), 0.0, Texture("data/mercury.jpg", false), nothing)
    push!(objs, Sphere(Vec3(-3, -2, -9) + offset, .383, mercury_material))

    venus_material = Material(Lambertian(), 0.0, Texture("data/venus.jpg", false), nothing)
    push!(objs, Sphere(Vec3(-1.2, -2, -8) + offset, .949, venus_material))
    
    earth_material = Material(Lambertian(), 0.0, Texture("data/better_earth.jpg", false), nothing)
    push!(objs, Sphere(Vec3(.5, -2, -6.5) + offset, 1, earth_material))
    
    mars_material = Material(Lambertian(), 0.0, Texture("data/mars.jpg", false), nothing)
    push!(objs, Sphere(Vec3(1.5, -2, -5) + offset, .532, mars_material))

    lights = [ DirectionalLight(0.8, Vec3(0,0,-1)),
               PointLight(0.8, Vec3(0,0,0)), 
            DirectionalLight(1, Vec3(-3, 2, 0))]

    air=[]
    Scene(bg, objs, lights, air)

end

#Final Project test scenes and helpers


function floor_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        temp = Vec3(mesh.positions[i][1]*(0.5), mesh.positions[i][2]*0.01, mesh.positions[i][3])
        mesh.positions[i] = temp * scale + translation
    end

    create_triangles(mesh, material)
end

function sidewall_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        temp = Vec3(mesh.positions[i][1]*0.01, mesh.positions[i][2], mesh.positions[i][3])
        mesh.positions[i] = temp * scale + translation
    end

    create_triangles(mesh, material)
end

function backwall_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        temp = Vec3(mesh.positions[i][1], mesh.positions[i][2], mesh.positions[i][3]*0.01)
        mesh.positions[i] = temp * scale + translation
    end

    create_triangles(mesh, material)
end

function windowtb_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        temp = Vec3(mesh.positions[i][1]*0.01, mesh.positions[i][2]*0.4, mesh.positions[i][3])
        mesh.positions[i] = temp * scale + translation
    end

    create_triangles(mesh, material)
end

function windows_helper(mesh, material, scale=1.0, translation=Vec3(0,0,0))

    for i in 1:length(mesh.positions)
        temp = Vec3(mesh.positions[i][1]*0.01, mesh.positions[i][2] * (0.5), mesh.positions[i][3]*0.4)
        mesh.positions[i] = temp * scale + translation
    end

    create_triangles(mesh, material)
end

function scene_12(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20)))
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20)))
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20)))
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -25)))

    #Window wall
    # append!(objs, sidewall_helper(cube_mesh(), cube_matg, 5, Vec3(5, 0, -20)))
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, -5, -20))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, 5, -20))) #top plane
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    # push!(objs, Sphere(Vec3(-2, -2, -19), 1.5, Material(Lambertian(), 0.2, nothing, blue))) 


    lights = [PointLight(white*2, Vec3(15, 3, -18))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(white*0.2, Vec3(0,0,-5))]
    # lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(white*0.2, Vec3(0,0,-5))]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), 1.5, .3)) 
    end   

    Scene(bg, objs, lights, air)
end

function scene_13(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    glossy = Material(Lambertian(), 0.0, 0.8, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20)))
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20)))
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20)))
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -25)))

    #Window wall
    # append!(objs, sidewall_helper(cube_mesh(), cube_matg, 5, Vec3(5, 0, -20)))
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, -5, -20))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, 5, -20))) #top plane
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    # push!(objs, Sphere(Vec3(-2, -2, -19), 1.5, Material(Lambertian(), 0.2, nothing, blue))) 


    lights = [PointLight(red*2, Vec3(15, 3, -18))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(white*0.2, Vec3(0,0,-5))]
    # lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(white*0.2, Vec3(0,0,-5))]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), 1.5, .3)) 
    end   

    Scene(bg, objs, lights, air)
end

function scene_14(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    glossy = Material(Lambertian(), 0.0, 0.8, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20)))
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20)))
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20)))
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -25)))

    #Window wall
    # append!(objs, sidewall_helper(cube_mesh(), cube_matg, 5, Vec3(5, 0, -20)))
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, -5, -20))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, 5, -20))) #top plane
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    # push!(objs, Sphere(Vec3(-2, -2, -19), 1.5, Material(Lambertian(), 0.2, nothing, blue))) 


    lights = [PointLight(green*2, Vec3(15, 3, -18))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(white*0.2, Vec3(0,0,-5))]
    # lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(white*0.2, Vec3(0,0,-5))]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), 1.5, .3)) 
    end   

    Scene(bg, objs, lights, air)
end

function scene_15(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    glossy = Material(Lambertian(), 0.0, 0.8, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20)))
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20)))
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20)))
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -25)))

    #Window wall
    # append!(objs, sidewall_helper(cube_mesh(), cube_matg, 5, Vec3(5, 0, -20)))
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, -5, -20))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, 5, -20))) #top plane
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    # push!(objs, Sphere(Vec3(-2, -2, -19), 1.5, Material(Lambertian(), 0.2, nothing, blue))) 


    lights = [PointLight(blue*2, Vec3(15, 3, -18))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(white*0.2, Vec3(0,0,-5))]
    # lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(white*0.2, Vec3(0,0,-5))]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), 1.5, .3)) 
    end   

    Scene(bg, objs, lights, air)
end

function scene_16(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    cube_matp = Material(Lambertian(), 0.0, nothing, purple)
    cube_matbk = Material(Lambertian(), 0.0, nothing, black)

    glossy = Material(Lambertian(), 0.5, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20))) #floor
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20))) #ceiling
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20))) #left wall
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -40)))

    #Window wall
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 14, Vec3(7, -8, -27))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 14, Vec3(7, 8, -27))) #top plane

    #segments to make windows - keep z values between -12 and -40
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -37))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #middle middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays in front of camera
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    push!(objs, Sphere(Vec3(-2, -2, -19), 0.7, Material(Lambertian(), 0.2, nothing, purple))) 
    cube_mat = Material(Lambertian(), 0.0, nothing, green)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(0, -6, -17)))
    cube_mat = Material(Lambertian(), 0.0, nothing, blue)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(-2, -6, -27)))
    cube_mat = Material(Lambertian(), 0.0, nothing, red)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(4, -6, -27)))

    lights = [PointLight(white, Vec3(15, 3, -18))]
    # lights = [PointLight(white, Vec3(15, 3, -18)), PointLight(white * .2, Vec3(0, 3, -5))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(red*0.2, Vec3(0,3,-5))]
    #lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(1.0, Vec3(0,0,-5)) ]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), .5, .3))
        push!(air, Cloud(Vec3(-2, -4, -22), 1.5, .5))  
    end   

    Scene(bg, objs, lights, air)
end

function scene_17(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    glossy = Material(Lambertian(), 0.0, 0.8, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20)))
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20)))
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20)))
    append!(objs, backwall_helper(cube_mesh(), glossy, 7, Vec3(0, 0, -25)))

    #Window wall
    # append!(objs, sidewall_helper(cube_mesh(), cube_matg, 5, Vec3(5, 0, -20)))
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, -5, -20))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, 5, -20))) #top plane
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    # push!(objs, Sphere(Vec3(-2, -2, -19), 1.5, Material(Lambertian(), 0.2, nothing, blue))) 


    lights = [PointLight(white*2, Vec3(15, 3, -18))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(white*0.2, Vec3(0,0,-5))]
    # lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(white*0.2, Vec3(0,0,-5))]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), 1.5, .3)) 
    end   

    Scene(bg, objs, lights, air)
end

function scene_18(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    cube_matp = Material(Lambertian(), 0.0, nothing, purple)
    cube_matbk = Material(Lambertian(), 0.0, nothing, black)
    glossy = Material(Lambertian(), 0.0, 0.3, nothing, white)


    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20))) #floor
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20))) #ceiling
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20))) #left wall
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -40)))

    #Window wall
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 14, Vec3(7, -8, -27))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 14, Vec3(7, 8, -27))) #top plane

    #segments to make windows - keep z values between -12 and -40
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -37))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #middle middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays in front of camera
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    push!(objs, Sphere(Vec3(-2, -2, -19), 0.7, Material(Lambertian(), 0.2, nothing, purple))) 
    append!(objs, mesh_helper(cube_mesh(), glossy, 0.8, Vec3(0, -6, -17)))
    cube_mat = Material(Lambertian(), 0.0, nothing, blue)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(-2, -6, -27))) 
    cube_mat = Material(Lambertian(), 0.0, nothing, red)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(4, -6, -27)))

    lights = [PointLight(white, Vec3(15, 3, -18))]
    #lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(1.0, Vec3(0,0,-5)) ]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), .5, .3))
        push!(air, Cloud(Vec3(-2, -4, -22), 1.5, .5))  
    end   

    Scene(bg, objs, lights, air)
end

function scene_19(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    glossy = Material(Lambertian(), 0.0, 0.3, nothing, white)
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, -7, -20)))
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20)))
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20)))
    append!(objs, backwall_helper(cube_mesh(), glossy, 7, Vec3(0, 0, -25)))

    #Window wall
    # append!(objs, sidewall_helper(cube_mesh(), cube_matg, 5, Vec3(5, 0, -20)))
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, -5, -20))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 7, Vec3(7, 5, -20))) #top plane
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    # push!(objs, Sphere(Vec3(-2, -2, -19), 1.5, Material(Lambertian(), 0.2, nothing, blue))) 


    lights = [PointLight(white*2, Vec3(15, 3, -18))]
    # lights = [PointLight(green, Vec3(15, 3, -18)), PointLight(white*0.2, Vec3(0,0,-5))]
    # lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(white*0.2, Vec3(0,0,-5))]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), 1.5, .3)) 
    end   

    Scene(bg, objs, lights, air)
end

function scene_20(local_fog::Bool)
    bg = black

    objs = []
    lights = []
    air = []

    cube_matb = Material(Lambertian(), 0.0, nothing, blue)
    cube_matr = Material(Lambertian(), 0.0, nothing, red)
    cube_matg = Material(Lambertian(), 0.0, nothing, green)
    cube_matw = Material(Lambertian(), 0.0, nothing, white)
    cube_matp = Material(Lambertian(), 0.0, nothing, purple)
    cube_matbk = Material(Lambertian(), 0.0, nothing, black)
    glossy = Material(Lambertian(), 0.0, 0.8, nothing, white)

    append!(objs, floor_helper(cube_mesh(), glossy, 20, Vec3(0, -7, -20))) #floor
    append!(objs, floor_helper(cube_mesh(), cube_matw, 20, Vec3(0, 7, -20))) #ceiling
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(-7, 0, -20))) #left wall
    append!(objs, backwall_helper(cube_mesh(), cube_matw, 7, Vec3(0, 0, -40)))

    #Window wall
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 14, Vec3(7, -8, -27))) #bottom plane
    append!(objs, windowtb_helper(cube_mesh(), cube_matg, 14, Vec3(7, 8, -27))) #top plane

    #segments to make windows - keep z values between -12 and -40
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -37))) #back middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -23))) #middle middle segment
    append!(objs, windows_helper(cube_mesh(), cube_matg, 7, Vec3(7, 0, -14))) #front middle segment

    #close off extra rays in front of camera
    append!(objs, sidewall_helper(cube_mesh(), cube_matr, 20, Vec3(8, 0, 5)))

    #room objects
    push!(objs, Sphere(Vec3(-5, -5, -19), 0.7, Material(Lambertian(), 0.0, nothing, blue))) 
    push!(objs, Sphere(Vec3(-2, -2, -19), 0.7, Material(Lambertian(), 0.2, nothing, purple))) 
    append!(objs, mesh_helper(cube_mesh(), glossy, 0.8, Vec3(0, -6, -17)))
    cube_mat = Material(Lambertian(), 0.0, nothing, blue)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(-2, -6, -27))) 
    cube_mat = Material(Lambertian(), 0.0, nothing, red)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.8, Vec3(4, -6, -27)))

    lights = [PointLight(white, Vec3(15, 3, -18))]
    #lights = [DirectionalLight(1.0, Vec3(1, 0.5, -0.1)), PointLight(1.0, Vec3(0,0,-5)) ]
   
    #preset localized density
    if local_fog
        push!(air, Cloud(Vec3(-2, -2, -19), .5, .3))
        push!(air, Cloud(Vec3(-2, -4, -22), 1.5, .5))  
    end   

    Scene(bg, objs, lights, air)
end

function artifact_laruek2(img_height, img_width)
    return Vec2(10,1)
end

function artifact_reada2(img_height, img_width)
    return [11, 1]
end

scenes = [scene_1, scene_2, scene_3, scene_4, scene_5, scene_6, scene_7, scene_8, scene_9, 
            scene_10, scene_11, scene_12, scene_13, scene_14, scene_15, scene_16, scene_17, scene_18, scene_19, scene_20]

end # module TestScenes
