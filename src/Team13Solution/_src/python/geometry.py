#!/usr/bin/python3

"""
    simulation of the Team Problem 13 - 3-D Non-Linear Magnetostatic Model
    
    12. Dez 2018 
    TU Wien
    
    valentin.hanser @ student.tuwien.ac.at
"""


def generateGeometry(air_box_size = [0.5, 0.5, 0.5], maxh_iron=1000, fullProblem=True):
    geo = CSGeometry()

    # set parameters
    core_width_outer= [0.2, 0.2]
    core_width_inner= [0.15, 0.15]
    core_height = 0.100
    
    height_sheets = 0.1264
    thickness_sheets = 0.0032




    # define cutting planes
    pl_top = Plane( Pnt(0, 0, core_height/2), Vec(0, 0, 1))
    pl_bottom = Plane( Pnt(0, 0, -core_height/2), Vec(0, 0, -1))

    # ------------------------------ 4 corners --------------------------------------
    # calc centres
    tmp_x = (core_width_inner[0]/2 - 0.025)
    tmp_y = (core_width_inner[1]/2 - 0.025)
    xy_centre = np.array([[tmp_x, tmp_y], [-tmp_x, tmp_y], [-tmp_x, -tmp_y], [tmp_x, -tmp_y]])

    # limiting planes
    pl_right_corner = Plane( Pnt(tmp_x, 0, 0), Vec(-1, 0, 0))
    pl_left_corner = Plane( Pnt(-tmp_x, 0, 0), Vec(1, 0, 0))
    pl_back_corner = Plane( Pnt(0, tmp_y, 0), Vec(0, -1, 0))
    pl_front_corner = Plane( Pnt(0, -tmp_y, 0), Vec(0, 1, 0))

    obj_stack = []

    # build 4 curves
    for i in range(4):
        if i == 0:
            # right back
            pl_a = pl_right_corner
            pl_b = pl_back_corner
            mat = "corner_right_back"
        elif i == 1:
            # left back
            pl_a = pl_left_corner
            pl_b = pl_back_corner
            mat = "corner_left_back"
        elif i == 2:
            # left front
            pl_a = pl_left_corner
            pl_b = pl_front_corner
            mat = "corner_left_front"
        elif i == 3:
            # right front
            pl_a = pl_right_corner
            pl_b = pl_front_corner
            mat = "corner_right_front"
        
        # create cylinders    
        cy_outer = Cylinder(Pnt(xy_centre[i, 0], xy_centre[i, 1], 0), Pnt(xy_centre[i,0], xy_centre[i, 1], 1), 0.05)
        cy_inner = Cylinder(Pnt(xy_centre[i, 0], xy_centre[i, 1], 0), Pnt(xy_centre[i, 0], xy_centre[i, 1], 1), 0.025)
         
        # create curved corner
        cy_ring = (cy_outer - cy_inner) * pl_top * pl_bottom * pl_a * pl_b
        # set properties
        cy_ring.mat(mat)
        # add to geo 
        obj_stack.append(cy_ring)
        
    # ------------------------------ bricks --------------------------------------
    pl_right_outer = Plane( Pnt(core_width_outer[0]/2, 0, 0), Vec(1, 0, 0))
    pl_left_outer = Plane( Pnt(-core_width_outer[0]/2, 0, 0), Vec(-1, 0, 0))
    pl_back_outer = Plane( Pnt(0, core_width_outer[1]/2, 0), Vec(0, 1, 0))
    pl_front_outer = Plane( Pnt(0, -core_width_outer[1]/2, 0), Vec(0, -1, 0))

    pl_right_inner = Plane( Pnt(core_width_inner[0]/2, 0, 0), Vec(-1, 0, 0))
    pl_left_inner = Plane( Pnt(-core_width_inner[0]/2, 0, 0), Vec(1, 0, 0))
    pl_back_inner = Plane( Pnt(0, core_width_inner[1]/2, 0), Vec(0, -1, 0))
    pl_front_inner = Plane( Pnt(0, -core_width_inner[1]/2, 0), Vec(0, 1, 0))

    # back
    brick_back = pl_top * pl_bottom * pl_back_outer * pl_back_inner - pl_right_corner - pl_left_corner
    brick_back.mat("brick_back")
    obj_stack.append(brick_back)

    # front
    brick_front = pl_top * pl_bottom * pl_front_outer * pl_front_inner - pl_right_corner - pl_left_corner
    brick_front.mat("brick_front")
    obj_stack.append(brick_front)

    # left
    brick_left = pl_top * pl_bottom * pl_left_outer * pl_left_inner - pl_back_corner - pl_front_corner
    brick_left.mat("brick_left")
    obj_stack.append(brick_left)

    # right
    brick_right = pl_top * pl_bottom * pl_right_outer * pl_right_inner - pl_back_corner - pl_front_corner
    brick_right.mat("brick_right")
    obj_stack.append(brick_right)

    # ------------------------------ vertical sheet --------------------------------------
    sheet_vertical = OrthoBrick (Pnt(-thickness_sheets/2, -0.05/2, -height_sheets/2), Pnt(thickness_sheets/2, 0.05/2, height_sheets/2))
    


    # ------------------------------ C sheets left --------------------------------------
    pl_C_left_back = Plane( Pnt(0, -0.04+0.05/2, 0), Vec(0, 1, 0))      # back plane
    pl_C_left_front = Plane( Pnt(0, -0.04-0.05/2, 0), Vec(0, -1, 0))    # front plane
    pl_C_left_right = Plane( Pnt(-0.0042/2, 0, 0), Vec(1, 0, 0))        # right plane

    sheet_C_left_outer = OrthoBrick (Pnt(-0.0042/2-0.12-thickness_sheets, -1, -height_sheets/2), Pnt(0, 1, height_sheets/2))              # outer block
    sheet_C_left_inner = OrthoBrick (Pnt(-0.0042/2-0.12, -1, -height_sheets/2 + thickness_sheets), Pnt(0, 1, height_sheets/2- thickness_sheets))    # inner block

    sheet_C_left = (sheet_C_left_outer - sheet_C_left_inner) * pl_C_left_back * pl_C_left_front * pl_C_left_right   # cut inner out of outer and cut around by planes


    # ------------------------------ C sheets right --------------------------------------
    pl_C_right_back = Plane( Pnt(0, 0.04+0.05/2, 0), Vec(0, 1, 0))
    pl_C_right_front = Plane( Pnt(0, 0.04-0.05/2, 0), Vec(0, -1, 0))
    pl_C_right_left = Plane( Pnt(0.0042/2, 0, 0), Vec(-1, 0, 0))

    sheet_C_right_outer = OrthoBrick (Pnt(0, -1, -height_sheets/2), Pnt(0.0042/2 + 0.12 + thickness_sheets, 1, height_sheets/2))
    sheet_C_right_inner = OrthoBrick (Pnt(0, -1, -height_sheets/2 + thickness_sheets), Pnt(0.0042/2 + 0.12, 1, height_sheets/2- thickness_sheets))

    sheet_C_right = (sheet_C_right_outer - sheet_C_right_inner) * pl_C_right_back * pl_C_right_front * pl_C_right_left


    sheets_iron = sheet_vertical + sheet_C_right + sheet_C_left
    sheets_iron.mat("iron")
    sheets_iron.maxh(maxh_iron)
    obj_stack.append(sheets_iron)
    
    
    # ----------------------- air box symmetries -----------------------------------
    if not fullProblem:
        pl_bottom = Plane( Pnt(0, 0, 0), Vec(0, 0,-1)).bc("mirror_z")
        
    else:
        pl_bottom = Plane( Pnt(0, 0, -air_box_size[2]/2), Vec(0, 0,-1)).bc("outer")

   
 
    pl_top = Plane(Pnt(0, 0, air_box_size[2]/2), Vec(0, 0, 1)).bc("outer")
    pl_left = Plane(Pnt(-air_box_size[0]/2, 0, 0), Vec(-1, 0, 0)).bc("outer")
    pl_right = Plane(Pnt(air_box_size[0]/2, 0, 0), Vec(1, 0, 0)).bc("outer")
    pl_back = Plane(Pnt(0, air_box_size[1]/2, 0), Vec(0, 1, 0)).bc("outer")
    pl_front = Plane(Pnt(0, -air_box_size[1]/2, 0), Vec(0, -1, 0)).bc("outer")
    
    air_box = pl_top * pl_left * pl_right * pl_back * pl_front * pl_bottom
        
    # subtract all items
    for i in range(len(obj_stack)):
        air_box -= obj_stack[i]
    air_box.mat("air")
   

    # --------------------------- add all items----------------------------
    geo.Add(air_box, maxh=100) #here for volume maxh
    for i in range(len(obj_stack)):
        geo.Add(obj_stack[i] * pl_bottom)
        
        
    

        
    return geo
    
def __restrictAreaZ(mp, x_range, y_range, z_range, Nx, Ny, Nz, maxh):
    xi = np.linspace(x_range[0], x_range[1], Nx)
    yi = np.linspace(y_range[0], y_range[1], Ny)
    zi = np.linspace(z_range[0], z_range[1], Nz)
    
    for z in zi:
        for y in yi:
            for x in xi:
                mp.RestrictH(x, y, z, h=maxh)
            
        maxh *=1.2
                
    return mp
    
def __restrictArea(mp, x_range, y_range, z_range, Nx, Ny, Nz, maxh, red="none"):
    xi = np.linspace(x_range[0], x_range[1], Nx)
    yi = np.linspace(y_range[0], y_range[1], Ny)
    zi = np.linspace(z_range[0], z_range[1], Nz)
    
    if red == "none":
        for x in xi:
            for y in yi:
                for z in zi:
                    mp.RestrictH(x, y, z, h=maxh)

                    
    elif red == "x":
        for x in xi:
            for y in yi:
                for z in zi:
                    mp.RestrictH(x, y, z, h=maxh)
                
            maxh *=1.2  
    elif red == "y":
        for y in yi:
            for x in xi:
                for z in zi:
                    mp.RestrictH(x, y, z, h=maxh)
                
            maxh *=1.2      
    elif red == "z":
        for z in zi:
            for x in xi:
                for y in yi:
                    mp.RestrictH(x, y, z, h=maxh)
                
            maxh *=1.2    
    return mp
    
    
def getMeshingParameters(maxh_global=100):
    mp = MeshingParameters(maxh=maxh_global)
    

    maxh_gap = 0.0032 / 4
    
    
    
    # vertical iron first area
    red = "z"
    mp = __restrictArea(mp, [0, 0], [0.01, 0.025], [0.126/2- 0.0032/2 , 0.126/2 - 0.03], 1, 3, 15, maxh_gap, red=red)
    mp = __restrictArea(mp, [0, 0], [-0.01, -0.025], [0.126/2- 0.0032/2 , 0.126/2 - 0.03], 1, 3, 15, maxh_gap, red=red)
    mp = __restrictArea(mp, [0, 0], [0.01, 0.025], [-0.126/2 + 0.0032/2 , -0.126/2 + 0.03], 1, 3, 15, maxh_gap, red=red)
    mp = __restrictArea(mp, [0, 0], [-0.01, -0.025], [-0.126/2 + 0.0032/2, -0.126/2 + 0.03], 1, 3, 15, maxh_gap, red=red)

    red = "x"
    # horizontal iron first area
    mp = __restrictArea(mp, [-0.0042/2, -0.03], [-0.015, -0.030], [0.126/2 - 0.0032/2 , 0.126/2 - 0.0032/2], 15, 3, 1, maxh_gap, red=red)
    mp = __restrictArea(mp, [-0.0042/2, -0.03], [-0.015, -0.030], [-0.126/2 + 0.0032/2 , -0.126/2 + 0.0032/2], 15, 3, 1, maxh_gap, red=red)
    mp = __restrictArea(mp, [0.0042/2, 0.03], [0.015, 0.030], [0.126/2 - 0.0032/2 , 0.126/2 - 0.0032/2], 15, 3, 1, maxh_gap, red=red)
    mp = __restrictArea(mp, [0.0042/2, 0.03], [0.015, 0.030], [-0.126/2 + 0.0032/2 , -0.126/2 + 0.0032/2], 15, 3, 1, maxh_gap, red=red)
    

    maxh_gap = 0.0032/1.5
    red = "x"
    # horizontal iron first area
    mp = __restrictArea(mp, [-0.0042/2, -0.07], [-0.04, -0.06], [0.126/2 - 0.0032 , 0.126/2], 7, 3, 3, maxh_gap, red=red)
    mp = __restrictArea(mp, [-0.0042/2, -0.07], [-0.04, -0.06], [-0.126/2  , -0.126/2 + 0.0032], 7, 3, 3, maxh_gap, red=red)
    mp = __restrictArea(mp, [0.0042/2, 0.07], [0.04, 0.06], [0.126/2 - 0.0032 , 0.126/2], 7, 3, 3, maxh_gap, red=red)
    mp = __restrictArea(mp, [0.0042/2, 0.07], [0.04, 0.06], [-0.126/2 + 0.0032 , -0.126/2], 7, 3, 3, maxh_gap, red=red)
    
    return mp




if __name__ == "__main__":
    from ngsolve import *
    from netgen.csg import *
    from netgen.meshing import MeshingParameters
    import sys, getopt
    import netgen.gui
    import numpy as np

    #default value
    fullProblem = False

    # input argument
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hf:",["fullProblem="])
    except getopt.GetoptError:
        print('geometry.py -f <bool fullProblem>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('geometry.py -f <bool fullProblem>')
            exit()
        elif opt in ("-f", "--fullProblem"):
            fullProblem = bool(arg)
    
    # generate geometry
    geo = generateGeometry(fullProblem=fullProblem)
    
    # meshing parameters    
    mp = getMeshingParameters()
    ngmesh = geo.GenerateMesh(mp=mp) # global maxh

    # save file
    ngmesh.Save("./team13_mesh.vol")

    # mesh it 
    mesh = Mesh(ngmesh)
    mesh.Curve(5)

    # draw domains
    val = {"corner_right_back":7.5 , "corner_left_back":2, "corner_left_front":3, "corner_right_front":4, "brick_front":6, "brick_back":5, "brick_left":7, "brick_right":8, "iron":9, "air":0}
    domains = CoefficientFunction([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])

    Draw(domains, mesh, "domains", draw_surf=False)
    input("finished")

else:

    from ngsolve import *
    from netgen.csg import *
    import numpy as np
    from netgen.meshing import MeshingParameters









