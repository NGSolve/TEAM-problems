from ngsolve import *
from netgen.csg import *

from myPackage import cmdInput, loadView, addExec
from geometry import generateGeometry, getMeshingParameters
from measurement import measure1to25, measure26to36, measure37to40
import os

import pickle
import datetime as dt

SetNumThreads(10)
import netgen.gui
import numpy as np
import matplotlib.pyplot as plt
import time

mu0 = 4*np.pi*1e-7
# ngsglobals.msg_level = 2
"""
    simulation of the Team Problem 13 - 3-D Non-Linear Magnetostatic Model
    
    12. Dez 2018 
    TU Wien
    
    valentin.hanser @ student.tuwien.ac.at
"""

# ------------------------- set parameters here -------------------------------
space_order = 2
I0N = 1000                  #Ampereturns

# ---------------------------- Meshing ----------------------------------------
# generate mesh
if False:
    geo = generateGeometry(maxh_iron = 100, fullProblem=False)
    mp = getMeshingParameters(maxh_global = 100)
    ngmesh = geo.GenerateMesh(mp=mp) # global maxh
    mesh = Mesh(ngmesh)
    mesh.Curve(5)
else:
    mesh = Mesh("./team13_mesh.vol")
    mesh.Curve(5)

# check domains
val = {"corner_right_back":7.5 , "corner_left_back":2, "corner_left_front":3, "corner_right_front":4, "brick_front":6, "brick_back":5, "brick_left":7, "brick_right":8, "iron":9, "air":0}
domains = CoefficientFunction([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
Draw(domains, mesh, "domains", draw_surf=False)

# --------------------------------- FE Space ----------------------------------
# create fe space
fes = HCurl(mesh, order=space_order, dirichlet="outer", nograds=True)

# magnetic vector potential as GridFunction
sol = GridFunction(fes, "A")
sol.vec[:] = 0

# Test- and Trialfunction
u = fes.TrialFunction()
v = fes.TestFunction()

# magnetic flux density
B = curl(sol)

B_GF = GridFunction(fes, "B")
B_GF.Set(B)
Draw(B, mesh, "B", draw_surf=False)

# ---------------------------------- material ---------------------------------
# calculate B/H curve according to paper

#H_KL, B_KL = createKL()
H_KL = [ -4.47197834e-13, 1.60000000e+01, 3.00000000e+01, 5.40000000e+01\
    , 9.30000000e+01, 1.43000000e+02, 1.91000000e+02, 2.10000000e+02\
, 2.22000000e+02, 2.33000000e+02, 2.47000000e+02, 2.58000000e+02\
, 2.72000000e+02, 2.89000000e+02, 3.13000000e+02, 3.42000000e+02\
, 3.77000000e+02, 4.33000000e+02, 5.09000000e+02, 6.48000000e+02\
, 9.33000000e+02, 1.22800000e+03, 1.93400000e+03, 2.91300000e+03\
, 4.99300000e+03, 7.18900000e+03, 9.42300000e+03, 9.42300000e+03\
, 1.28203768e+04, 1.65447489e+04, 2.07163957e+04, 2.55500961e+04\
, 3.15206135e+04, 4.03204637e+04, 7.73038295e+04, 1.29272791e+05\
, 1.81241752e+05, 2.33210713e+05, 2.85179674e+05, 3.37148635e+05\
, 3.89117596e+05, 4.41086557e+05, 4.93055518e+05, 5.45024479e+05\
, 5.96993440e+05, 6.48962401e+05, 7.00931362e+05, 7.52900323e+05\
, 8.04869284e+05, 8.56838245e+05, 9.08807206e+05, 9.60776167e+05\
, 1.01274513e+06, 1.06471409e+06, 1.11668305e+06, 1.16865201e+06\
, 1.22062097e+06, 1.27258993e+06, 1.32455889e+06, 1.37652785e+06\
, 1.42849682e+06, 1.48046578e+06, 1.53243474e+06, 1.58440370e+06\
, 1.63637266e+06, 1.68834162e+06, 1.74031058e+06, 1.79227954e+06\
, 1.84424850e+06, 1.89621746e+06, 1.94818643e+06, 2.00015539e+06\
, 2.05212435e+06, 2.10409331e+06, 2.15606227e+06, 2.20803123e+06\
, 2.26000019e+06]

B_KL = [  0.00000000e+00, 2.50000000e-03, 5.00000000e-03, 1.25000000e-02\
, 2.50000000e-02, 5.00000000e-02, 1.00000000e-01, 2.00000000e-01\
, 3.00000000e-01, 4.00000000e-01, 5.00000000e-01, 6.00000000e-01\
, 7.00000000e-01, 8.00000000e-01, 9.00000000e-01, 1.00000000e+00\
, 1.10000000e+00, 1.20000000e+00, 1.30000000e+00, 1.40000000e+00\
, 1.50000000e+00, 1.55000000e+00, 1.60000000e+00, 1.65000000e+00\
, 1.70000000e+00, 1.75000000e+00, 1.80000000e+00, 1.80000000e+00\
, 1.86530612e+00, 1.93061224e+00, 1.99591837e+00, 2.06122449e+00\
, 2.12653061e+00, 2.19183673e+00, 2.25714286e+00, 2.32244898e+00\
, 2.38775510e+00, 2.45306122e+00, 2.51836735e+00, 2.58367347e+00\
, 2.64897959e+00, 2.71428571e+00, 2.77959184e+00, 2.84489796e+00\
, 2.91020408e+00, 2.97551020e+00, 3.04081633e+00, 3.10612245e+00\
, 3.17142857e+00, 3.23673469e+00, 3.30204082e+00, 3.36734694e+00\
, 3.43265306e+00, 3.49795918e+00, 3.56326531e+00, 3.62857143e+00\
, 3.69387755e+00, 3.75918367e+00, 3.82448980e+00, 3.88979592e+00\
, 3.95510204e+00, 4.02040816e+00, 4.08571429e+00, 4.15102041e+00\
, 4.21632653e+00, 4.28163265e+00, 4.34693878e+00, 4.41224490e+00\
, 4.47755102e+00, 4.54285714e+00, 4.60816327e+00, 4.67346939e+00\
, 4.73877551e+00, 4.80408163e+00, 4.86938776e+00, 4.93469388e+00\
, 5.00000000e+00]
bh_curve = BSpline (2, [0]+list(B_KL), list(H_KL))

energy = bh_curve.Integrate()    # to minimise 

# ----------------------------------- lhs -------------------------------------
a = BilinearForm(fes, symmetric=False)

a += SymbolicBFI(1/mu0 * curl(u)*curl(v), definedon=~mesh.Materials("iron"))
# a += SymbolicBFI(1/(mu0*mur_iron) * curl(u)*curl(v), definedon=mesh.Materials("iron")) # linear
a += SymbolicEnergy(energy(sqrt(1e-12+curl(u)*curl(u))), definedon=mesh.Materials("iron")) # 1e-12 ... regularisation, avoid 0
a += SymbolicBFI(1e-1*u*v)

c = Preconditioner(a, type="direct", inverse="sparsecholesky")

# ------------------------------------ rhs ------------------------------------
A = 2500*1e-6       
J0 = I0N/A
# +++ bricks +++
J_brick_back = [1, 0]
J_brick_front = [-1, 0]
J_brick_left = [0, 1]
J_brick_right = [0, -1]

# +++ corners +++
# right back
x_right_back = x - 0.050
y_right_back = y - 0.050
r_right_back = (x_right_back**2 + y_right_back**2)**(1/2)
J_corner_right_back = [1/r_right_back * y_right_back, -1/r_right_back * x_right_back]

# left back
x_left_back = x + 0.050
y_left_back = y - 0.050
r_left_back = (x_left_back**2 + y_left_back**2)**(1/2)
J_corner_left_back =  [1/r_left_back * y_left_back, -1/r_left_back * x_left_back]


# left front
x_left_front = x + 0.050
y_left_front = y + 0.050
r_left_front = (x_left_front**2 + y_left_front**2)**(1/2)
J_corner_left_front = [1/r_left_front * y_left_front, -1/r_left_front * x_left_front]

# right front
x_right_front = x - 0.050
y_right_front = y + 0.050
r_right_front = (x_right_front**2 + y_right_front**2)**(1/2)
J_corner_right_front = [1/r_right_front * y_right_front, -1/r_right_front * x_right_front]

# r
val={"corner_right_back":r_right_back, "corner_left_back":r_left_back, "corner_left_front":r_left_front, "corner_right_front":r_right_front}
radius = CoefficientFunction([val[mat] if mat in val.keys() else 0 for mat in mesh.GetMaterials()])
Draw(radius, mesh, "radius")

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

# -------------------------------------- solve --------------------------------
t_simulation = time.time()
with TaskManager():
    f.Assemble()
    B_GF.Update()
    
    err = 1
    it = 1

    au = sol.vec.CreateVector()
    r = sol.vec.CreateVector()
    w = sol.vec.CreateVector()
    sol_new = sol.vec.CreateVector()

    while err > 1e-10:
        print ("nonlinear iteration", it)
        it = it+1

        E0 = a.Energy(sol.vec) - InnerProduct(f.vec, sol.vec)
        print ("Energy old = ", E0)
        
        a.AssembleLinearization(sol.vec)
        a.Apply (sol.vec, au)
        r.data = f.vec - au

        inv = CGSolver (mat=a.mat, pre=c.mat)
        w.data = inv * r

        err = InnerProduct (w, r)
        print ("err = ", err)


        sol_new.data = sol.vec + w
        E = a.Energy(sol_new) - InnerProduct(f.vec, sol_new)
        print ("Enew = ", E)
        tau = 1
        while E > E0:
            tau = 0.5*tau
            sol_new.data = sol.vec + tau * w
            E = a.Energy(sol_new) - InnerProduct(f.vec, sol_new)
            print ("tau = ", tau, "Enew =", E)

        sol.vec.data = sol_new

        Redraw()
t_simulation = time.time() - t_simulation
print("simulation time %.3lf seconds" % t_simulation)

# ---------------------------------- measurements -----------------------------
# results according to Team Problem 13 description page 13
B_msm = [   1.33, 1.329, 1.286, 1.225, 1.129, 0.985, 0.655, \
            0.259, 0.453, 0.554, 0.637, 0.698, 0.755, 0.809, 0.901, 0.945, 0.954, 0.956,\
            0.960, 0.965, 0.970, 0.974, 0.981, 0.984, 0.985]
            
# measure avarage magnetic flux density            
B_sim_fir = measure1to25(B, mesh, draw=False)
B_sim_sec = measure26to36(B, mesh)
B_sim_thi = measure37to40(B, mesh)

B_sim = np.hstack([B_sim_fir, B_sim_sec, B_sim_thi])

# output
for i in range(len(B_msm)):
    print("measurement %d: sim\t %.3lf, msm\t %.3lf, rel. err\t %.3lf" % (i + 1, B_sim[i], B_msm[i], (B_sim[i] - B_msm[i])/ B_msm[i] * 100 ))

for i in range(i + 1, len(B_sim)):
    print("measurement %d: sim\t %.3lf" % (i + 1, B_sim[i]))
    
# -------------------------------- create a plot ------------------------------
# plot magnetic flux density over position 
plt.figure(1, figsize=[12,9])
plt.clf()
xi = range(0, 25)
plt.plot([7, 7], [0, 1.6], "--k")
plt.plot([18, 18], [0, 1.6], "--k")
plt.plot(xi, B_msm[0:25], '--o', label="measured")
plt.plot(xi, B_sim[0:25], '--D', label="simulated")
plt.title("comparison of measured and simulated values")
plt.ylabel("according component of B in T")
plt.xlabel("measurement number")
plt.ylim(0, 1.6)
plt.grid()

plt.annotate("air gap",xy=(7, 1), arrowprops=dict(arrowstyle='->'), xytext=(10, 1.2))
plt.annotate("corner",xy=(18, 0.2), arrowprops=dict(arrowstyle='->'), xytext=(15, 0.6))

plt.legend(loc=0)
font = {'size'   : 16}
matplotlib.rc('font', **font)
plt.show(block=False)

# ---------------------------------- finish -----------------------------------
# draw Bnorm
fesBnorm = H1(mesh, definedon=mesh.Materials("iron"))
B_norm = GridFunction(fesBnorm, "|B|")
B_norm.Set(B.Norm())
Draw(B_norm)

# save to a file
toSave={"fes":fes, "A":sol.vec, "B":B_GF.vec}
if not os.path.exists("./save"):
    os.mkdir("./save")
pickle.dump(toSave, open( "./save/teamProblem13_"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S") +"_order"+str(space_order)+"_ndofs_"+str(fes.ndof)+".pkl", "wb" ))

print("ndofs: %d s" % (fes.ndof))
print("running time of simulation: %.3lf s" % (t_simulation))
loadView()
cmdInput(locals(), globals())


