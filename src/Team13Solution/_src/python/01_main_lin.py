from myPackage import cmdInput, loadView
from ngsolve import *
from netgen.csg import *

from geometry import generateGeometry

from netgen.meshing import MeshingParameters

SetNumThreads(10)
import netgen.gui

import numpy as np

#------------------------------------------------

space_order = 1

mur_iron = 1e4
mu0 = 4*np.pi*1e-7

J0 = 1000 

#------------------------------------------------

# generate mesh
geo = generateGeometry(maxh_iron=0.01)


mp = MeshingParameters(maxh=100)
mp.RestrictH(-0.003, -0.015 + i, 0.1264/2, h=0.001)


ngmesh = geo.GenerateMesh(mp=mp) # global maxh
mesh = Mesh(ngmesh)
mesh.Curve(5)

val = {"corner_right_back":7.5 , "corner_left_back":2, "corner_left_front":3, "corner_right_front":4, "brick_front":6, "brick_back":5, "brick_left":7, "brick_right":8, "iron":9, "air":0}
domains = CoefficientFunction([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
Draw(domains, mesh, "domains", draw_surf=False)


val = {"outer":1 , "default":0}
boundaries = CoefficientFunction([val[bnd] if bnd in val.keys() else 0 for bnd in mesh.GetBoundaries()])
Draw(boundaries, mesh, "boundaries", draw_surf=True)

# ------------------------ solve --------------------------------
fes = HCurl(mesh, order=space_order, dirichlet="outer", nograds=True)

sol = GridFunction(fes, "A")

u = fes.TrialFunction()
v = fes.TestFunction()

B = curl(sol)


# ------------------------- lhs ----------------------------------
a = BilinearForm(fes, symmetric=False)

a += SymbolicBFI(1/mu0 * curl(u)*curl(v), definedon=~mesh.Materials("iron"))
a += SymbolicBFI(1/(mu0*mur_iron) * curl(u)*curl(v), definedon=mesh.Materials("iron"))
a += SymbolicBFI(1e-1*u*v)


# ------------------------- rhs ---------------------------------- 
J_brick_back = [1, 0]
J_brick_front = [-1, 0]
J_brick_left = [0, 1]
J_brick_right = [0, -1]

# right back
x_right_back = x - 0.050
y_right_back = y - 0.050
r_right_back = (x_right_back**2 + y_right_back**2)**(1/2)
# phi_right_back = acos(x_right_back/r_right_back)
J_corner_right_back = [1/r_right_back * y_right_back, -1/r_right_back * x_right_back]

# left back
x_left_back = x + 0.050
y_left_back = y - 0.050
r_left_back = (x_left_back**2 + y_left_back**2)**(1/2)
# phi_left_back = acos(x_left_back/r_left_back)
J_corner_left_back =  [1/r_left_back * y_left_back, -1/r_left_back * x_left_back]


# left front
x_left_front = x + 0.050
y_left_front = y + 0.050
r_left_front = (x_left_front**2 + y_left_front**2)**(1/2)
# phi_left_front = 2*np.pi - acos(x_left_front/r_left_front)
J_corner_left_front = [1/r_left_front * y_left_front, -1/r_left_front * x_left_front]



# right front
x_right_front = x - 0.050
y_right_front = y + 0.050
r_right_front = (x_right_front**2 + y_right_front**2)**(1/2)
# phi_right_front = 2*np.pi - acos(x_right_front/r_right_front)
J_corner_right_front = [1/r_right_front * y_right_front, -1/r_right_front * x_right_front]


# r
val={"corner_right_back":r_right_back, "corner_left_back":r_left_back, "corner_left_front":r_left_front, "corner_right_front":r_right_front}
radius = CoefficientFunction([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
Draw(radius, mesh, "radius")

# phi
# val={"corner_right_back":phi_right_back, "corner_left_back":phi_left_back, "corner_left_front":phi_left_front, "corner_right_front":phi_right_front}
# phi = CoefficientFunction([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
# Draw(phi, mesh, "phi")


# J

val={   "corner_right_back":J_corner_right_back, "corner_left_back":J_corner_left_back, \
        "corner_left_front":J_corner_left_front, "corner_right_front":J_corner_right_front,\
        "brick_back":J_brick_back, "brick_left":J_brick_left, \
        "brick_front":J_brick_front, "brick_right":J_brick_right}    


J = J0 * CoefficientFunction([val[mat][0] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])*CoefficientFunction((1,0,0)) + \
    J0 * CoefficientFunction([val[mat][1] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])*CoefficientFunction((0,1,0))

Draw(J, mesh, "J")


f = LinearForm(fes)
f += SymbolicLFI(J*v)

# -------------------------------------- solve -----------------------------
with TaskManager():
    c = Preconditioner(a, type = "direct")
    f.Assemble()
    a.Assemble()
    bvp = BVP(bf = a, lf = f, gf = sol, pre = c, maxsteps=50)
    bvp.Do()
    
Draw(B, mesh, "B", draw_surf=False)

loadView()
cmdInput(locals(), globals())


