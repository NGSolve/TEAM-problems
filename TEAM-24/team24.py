from geometry import MakeGeometry
from netgen.meshing import *
from ngsolve import *
from ngsolve.nonlinearsolvers import NewtonSolver
from ngsolve.krylovspace import CGSolver
from math import pi
from bhcurve import GetBHCurve
import matplotlib
matplotlib.use("WXAGG")
from hall_point import hall_point

ngsglobals.msg_level = 1

SetNumThreads(6)
SetHeapSize(10**8)

mp = MeshingParameters(meshsize.moderate)
mp.RestrictH(*hall_point, 0.0001)
mp.RestrictH(-hall_point[0], -hall_point[1], hall_point[2], 0.0001)

geo = MakeGeometry()
with TaskManager():
    m1 = geo.GenerateMesh(mp, maxh=0.02)
    mesh = Mesh(m1)
    mesh.Curve(5)

Draw(mesh)

hall_probe = mesh(*hall_point)

coilcur = { "coil1a" : (0,0,1),
            "coil2a" : (z-0.0155, 0, -x-0.017),
            "coil3a" : (1,0,0),
            "coil4a" : (z-0.0155, 0, 0.017-x),
            "coil5a" : (0,0,-1),
            "coil6a" : (z+0.0155, 0, 0.017-x),
            "coil7a" : (-1,0,0),
            "coil8a" : (z+0.0155, 0, -0.017-x),
            "coil1b" : (0,0,1),
            "coil2b" : (z-0.0155, 0, -x-0.017),
            "coil3b" : (1,0,0),
            "coil4b" : (z-0.0155, 0, 0.017-x),
            "coil5b" : (0,0,-1),
            "coil6b" : (z+0.0155, 0, 0.017-x),
            "coil7b" : (-1,0,0),
            "coil8b" : (z+0.0155, 0, -0.017-x)
        }
cur = [ coilcur[mat] if mat in coilcur else (0,0,0) for mat in mesh.GetMaterials() ]
cur = -CoefficientFunction(cur)
cur = CoefficientFunction([1.0/(Norm(cur)+1e-20) * cur if mat in coilcur else (0,0,0) for mat in mesh.GetMaterials()])

coildomains = mesh.Materials("coil[1-9][a-b]")
Draw (cur, mesh, "current")

graddoms = [ 1 if mat=="stator" or mat == "rotor" else 0 for mat in mesh.GetMaterials() ]
order = 3
V = HCurl(mesh, order=order, dirichlet="sym", gradientdomains=graddoms)
print ("coil - domains = ", coildomains.Mask())
VI = NumberSpace(mesh)
X = FESpace([V,VI])

u,I = X.TrialFunction()
v,Itilde = X.TestFunction()

mu0 = 1.257e-6
nu = 1./mu0

bhcurve = GetBHCurve()

magnetic = mesh.Materials("stator|rotor")
print("magnetic: ", magnetic.Mask())

uold = GridFunction(V)
tau = 0.005
tend = 0.3-1e-6
# tend = 3*tau

sigma = CoefficientFunction ([ 4.54e5 if mat=="stator" or mat == "rotor" else 0 for mat in mesh.GetMaterials() ])
0.0139, 0.0519216
def normeps(f):
    return sqrt(f*f + 1e-20)

R = 3.09 # Ohm

N = 350
S = 0.024*0.017  # coil cross-section
B = curl(u)
nu = CoefficientFunction([bhcurve(normeps(B))/normeps(B) if mat in ["stator", "rotor"] else 1./mu0 for mat in mesh.GetMaterials()])

a = BilinearForm(X)
a += (nu * (curl(u) * curl(v) + 1e-6 * u * v) + sigma/tau * (u-uold) * v) * dx
a += -N/S * I * cur * v * dx(coildomains, bonus_intorder=30)

voldom = Integrate(1, mesh, definedon=coildomains)
# actually \int E = \int (u-uold)/tau must be multiplied by 2
# multiplying the voldom is the same and keeps symmetry
voldom *= 2 # for symmetry

U = 23.1
a += tau/voldom * (-R*I-U) * Itilde * dx(coildomains)

a += - N/S * cur * (u-uold) * Itilde * dx(coildomains, bonus_intorder=30)
if order < 2:
    c = Preconditioner(a, "direct", inverse="sparsecholesky")
else:
    # c = Preconditioner(a,"bddc", coarsetype="myamg_hcurl")
    c = Preconditioner(a, "bddc", inverse="sparsecholesky")

u = GridFunction(X)
uprime = GridFunction(V)

SetVisualization(clipping=True, clipnormal=(0,0,-1))
Draw (curl(uold), mesh, "B-field")
Draw (sigma * uprime, mesh, "Eddy-current")
from ngsolve.internal import visoptions, viewoptions
visoptions.clipsolution = "vec"
viewoptions.clipping.dist = 0.190
visoptions.gridsize = 120

uold.vec[:] = 0.0

res = u.vec.CreateVector()
w = u.vec.CreateVector()

times = [0]
By_vals = [0]
I_vals = [0]
rotor_flux = [0]

measured_t = [0, 0.002, 0.005, 0.010, 0.015, 0.020,
              0.030, 0.040, 0.050, 0.060, 0.070, 0.080,
              0.090, 0.100, 0.120, 0.140, 0.160, 0.180,
              0.200, 0.240, 0.300]
measured_By = [0, 0.07, 0.21, 0.38, 0.53, 0.64, 0.82, 0.94, 1.03,
               1.09, 1.14, 1.17, 1.19, 1.21, 1.23, 1.24, 1.25,
               1.25, 1.25, 1.26, 1.26]
measured_I = [0, 0.91, 1.80, 2.95, 3.82, 4.49, 5.45, 6.05, 6.45,
              6.71, 6.91, 7.04, 7.15, 7.22, 7.31, 7.36, 7.38, 7.40,
              7.40, 7.41, 7.41]
measured_rotor_flux = [0, 0.29, 0.77, 1.44, 1.99, 2.44, 3.11, 3.58, 3.91, 4.14, 4.31, 4.43, 4.52, 4.58, 4.65, 4.69, 4.71, 4.72, 4.72, 4.72, 4.72]
measured_rotor_flux = [1e-4*val for val in measured_rotor_flux]

import matplotlib.pyplot as plt
nplots = 3
i = 0
plt.ion()
fig, axs = plt.subplots(nplots, 1, constrained_layout=True)
# relative error
plot_by, = axs[i].plot(times, By_vals, label="computed")
axs[i].plot(measured_t, measured_By, label="measured")
axs[i].set_title("By at hall probe")
axs[i].legend()
i += 1

plot_i, = axs[i].plot(times, I_vals, label="computed")
axs[i].plot(measured_t, measured_I, label="measured")
axs[i].set_title("I")
axs[i].legend()
i += 1

plot_flux, = axs[i].plot(times, rotor_flux, label="computed")
axs[i].plot(measured_t, measured_rotor_flux, label="measured")
axs[i].set_title("Rotor flux")
axs[i].legend()
i += 1

fig.canvas.draw()
fig.canvas.flush_events()

with TaskManager():
    t = 0
    a.AssembleLinearization(u.vec)
    solver = CGSolver(a.mat, c, maxsteps=1000, abstol=1e-10)
    newton = NewtonSolver(a, u, solver=solver)
    while t < tend:
        t = t + tau
        print("t = ", t)
        newton.Solve(maxerr=1e-8, printing=True)
        uprime.vec.data = 1/tau * u.components[0].vec - 1/tau * uold.vec
        uold.vec.data = u.components[0].vec
        times.append(t)
        by = curl(uold)(hall_probe)[1]
        I = -u.components[1].vec[0]
        print("By = ", by)
        print("I = ", I)

        n = specialcf.normal(3)
        BA = Integrate(curl(uold) * IfPos(y, 1, -1) * specialcf.normal(3), mesh, definedon=mesh.Boundaries("integrate"))/2
        B_rotor = Integrate(curl(uold) * IfPos(y, 1, 0) * n, mesh, definedon=mesh.Boundaries("rotor_end"))
        print("rotor flux = ", B_rotor)
        rotor_flux.append(B_rotor)
        By_vals.append(by)
        I_vals.append(I)
        for plot in (plot_by, plot_i, plot_flux):
            plot.set_xdata(times)
        plot_by.set_ydata(By_vals)
        plot_i.set_ydata(I_vals)
        plot_flux.set_ydata(rotor_flux)
        fig.canvas.draw()
        fig.canvas.flush_events()
        Redraw()


print("By = ", By_vals)


plt.savefig("results.png")

