from netgen.occ import *
from ngsolve import *
from ngsolve.krylovspace import CG

ngsglobals.msg_level = 0

configuration = "A"

if configuration == "A":
    d1 = 0.003048
    d2 = 0.0039624
    d3 = 0.0029972
    l = 0.0016
    lc = 0.001524
    sc = 0.000127
    delta = 0.000254
elif configuration == "B":
    d1 = 0.001524
    d2 = 0.003175
    d3 = 0.0016
    l = 0.000818
    lc = 0.001524    
    sc = 0.000127
    delta = 0.000254

r = 0.01

magnet = Cylinder((delta, 0, 0), X, d3/2, l)
magnet.faces.name = "bnd_magnet"
magnet.solids.name = "magnet"

coil = Cylinder((-sc, 0, 0), -X, d2/2, lc) - \
    Cylinder((-sc, 0, 0), -X, d1/2, lc)
coil.faces.name = "bnd_coil"
coil.solids.name = "coil"

box = Box((-r,-r,-r), (r,r,r))
air = box - coil - magnet
air.solids.name = "air"

shape = Glue([magnet, coil, air])

geo = OCCGeometry(shape)

mesh = Mesh(geo.GenerateMesh(maxh=0.01))
# mesh.Refine()
# mesh.Refine()
mesh.Curve(3)
Draw(mesh)

fes = HCurl(mesh, order=3, nograds=True)
u,v = fes.TnT()

mu0 = 1.257e-6
mur = mesh.MaterialCF({ "magnet" : 1.22 }, default=1)
nu = 1/(mu0*mur)

a = BilinearForm(fes)
a += nu * (curl(u) * curl(v) + 1e-6 * u * v) *dx
c = Preconditioner(a, "bddc")

f = LinearForm(fes)
tau = CoefficientFunction((0, z, -y))
tau *= 1/Norm(tau)
J = 0.1 * 280 / (lc * (d2-d1)/2) * tau # Amp turns / coil area * tau
f += J * v * dx("coil")

gfu = GridFunction(fes)

with TaskManager():
    a.Assemble()
    f.Assemble()
    solver = CGSolver(mat=a.mat, pre=c, printrates=True)
    gfu.vec.data = solver * f.vec

    B = curl(gfu)
    Draw(B, mesh, "B")

    # Maxwell stress tensor
    sigma = nu * (OuterProduct(B,B) - 0.5 * (B*B) * Id(3))

    # Compute fore in H1-space
    fes = VectorH1(mesh, order=3)
    u,v = fes.TnT()
    force = LinearForm(fes)
    force += InnerProduct(sigma, grad(v)) * dx
    force.Assemble()

    mass = BilinearForm(fes, check_unused=False)
    mass += u * v * ds
    mass.Assemble()

    h1force = GridFunction(fes)
    h1force.vec.data = mass.mat.Inverse(freedofs=fes.GetDofs(mesh.Boundaries("bnd_magnet|bnd_coil"))) * force.vec

    # compute surface force in Facet-space as  [sigma n]
    fes = FacetFESpace(mesh, order=3)**3
    u,v = fes.TnT()
    n = specialcf.normal(3)

    force = LinearForm(fes)
    force += (sigma * n) * v * dx(element_boundary=True)
    force.Assemble()

    mass = BilinearForm(fes, check_unused=False)
    mass += 0.5 * u * v * dx(element_boundary=True)
    mass.Assemble()

    facetforce = GridFunction(fes)
    facetforce.vec.data = mass.mat.Inverse(freedofs=fes.GetDofs(mesh.Boundaries("bnd_magnet|bnd_coil"))) * force.vec

    print("force-h1, magnet:", Integrate(h1force, mesh.Boundaries("bnd_magnet")))
    print("force-facet, magnet:", Integrate(facetforce, mesh.Boundaries("bnd_magnet")))

    print("force-h1, coil:", Integrate(h1force, mesh.Boundaries("bnd_coil")))
    print("force-facet, coil:", Integrate(facetforce, mesh.Boundaries("bnd_coil")))

    print("Lorentz force on coil:", Integrate(Cross(B,J), mesh.Materials("coil")))
