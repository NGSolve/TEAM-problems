from netgen.csg import *
from ngsolve import *
import numpy as np
from myPackage import cmdInput


def MakeGeometry():
    geometry = CSGeometry()
    box = OrthoBrick(Pnt(-1,-1,-1),Pnt(2,1,2)).bc("outer")

    core = OrthoBrick(Pnt(0,-0.05,0),Pnt(0.8,0.05,1))- \
           OrthoBrick(Pnt(0.1,-1,0.1),Pnt(0.7,1,0.9)).mat("core")
    core.maxh(0.05)
    
    coil = (Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.3) - \
            Cylinder(Pnt(0.05,0,0), Pnt(0.05,0,1), 0.15)) * \
            OrthoBrick (Pnt(-1,-1,0.3),Pnt(1,1,0.7)).maxh(0.2).mat("coil")

    geometry.Add ((box-core-coil).mat("air"))
    geometry.Add (core)
    geometry.Add (coil)
    return geometry




mesh = Mesh (MakeGeometry().GenerateMesh(maxh=100))#0.5))
mesh.Curve(5)

# ngsglobals.msg_level = 5


V = HCurl(mesh, order=1, dirichlet="outer", nograds=True)
u,v = V.TnT()

bhfilename = "bh-105-30h_50hz_TKES.dat"
f = open(bhfilename, 'r')

bhvalues = []
num = int(f.readline())
print ("num = ", num)
for line in f:
    b,h = line.split()
    bhvalues.append ( (float(b),float(h)) )

b,h = zip (*bhvalues)
print ("b=", b, "h = ", h)



bhcurve = BSpline (2, [0]+list(b), list(h))
energy = bhcurve.Integrate()


mu0 = 1.257e-6
mur = { "core" : 1000, "coil" : 1, "air" : 1 }
nu_coef = [ 1/(mu0*mur[mat]) for mat in mesh.GetMaterials() ]

nu = CoefficientFunction(nu_coef)
a = BilinearForm(V, symmetric=True)
a += SymbolicBFI(nu*curl(u)*curl(v), definedon=~mesh.Materials("core"))
a += SymbolicEnergy(energy(sqrt(1e-12+curl(u)*curl(u))), definedon=mesh.Materials("core"))
a += SymbolicBFI(1e-6/mu0*u*v)

c = Preconditioner(a, type="direct", inverse="sparsecholesky")


f = LinearForm(V)
f += SymbolicLFI(1e4*CoefficientFunction((y,-x,0)) * v, definedon=mesh.Materials("coil"))
f.Assemble()

u = GridFunction(V)
u.vec[:] = 0



Draw (curl(u), mesh, "B-field", draw_surf=False)



#  with TaskManager():
#    solvers.NewtonMinimization(a, u, linesearch=True)


err = 1
it = 1

au = u.vec.CreateVector()
r = u.vec.CreateVector()
w = u.vec.CreateVector()
unew = u.vec.CreateVector()

while err > 1e-10:
    print ("nonlinear iteration", it)
    it = it+1

    E0 = a.Energy(u.vec) - InnerProduct(f.vec, u.vec)
    print ("Energy old = ", E0)
    
    a.AssembleLinearization(u.vec)
    a.Apply (u.vec, au)
    r.data = f.vec - au

    inv = CGSolver (mat=a.mat, pre=c.mat)
    w.data = inv * r

    err = InnerProduct (w, r)
    print ("err = ", err)


    unew.data = u.vec + w
    E = a.Energy(unew) - InnerProduct(f.vec, unew)
    print ("Enew = ", E)
    tau = 1
    while E > E0:
        tau = 0.5*tau
        unew.data = u.vec + tau * w
        E = a.Energy(unew) - InnerProduct(f.vec, unew)
        print ("tau = ", tau, "Enew =", E)

    u.vec.data = unew

    Redraw()
    # input ("<waiting for you>")


