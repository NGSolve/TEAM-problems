from myPackage import cmdInput, loadView

from ngsolve import *
import netgen.gui
import pickle
import numpy as np

from measurement import measure1to25, measure26to36, measure37to40
import matplotlib.pyplot as plt

# load stuff
loaded = pickle.load(open( "./save/teamProblem13_2018-12-12 15:43:46_order2_ndofs_194296.pkl", "rb" ))

fes = loaded["fes"]
A = GridFunction(fes, "A")
A.vec.data = loaded["A"]


B = curl(A)
# B = GridFunction(fes, "B")
# B.vec.data = loaded["B"]


mesh = fes.mesh
# draw Bnorm
fesBnorm = H1(mesh, definedon=mesh.Materials("iron"))
B_norm = GridFunction(fesBnorm, "|B|")
B_norm.Set(B.Norm())
Draw(B_norm)


# print failure relative to measurement
B_msm = [   1.33, 1.329, 1.286, 1.225, 1.129, 0.985, 0.655, \
            0.259, 0.453, 0.554, 0.637, 0.698, 0.755, 0.809, 0.901, 0.945, 0.954, 0.956,\
            0.960, 0.965, 0.970, 0.974, 0.981, 0.984, 0.985]
B_sim_fir = measure1to25(B, mesh, draw=False)
B_sim_sec = measure26to36(B, mesh)
B_sim_thi = measure37to40(B, mesh)

B_sim = np.hstack([B_sim_fir, B_sim_sec, B_sim_thi])
# print results
for i in range(len(B_msm)):
    print("measurement %d: sim\t %.3lf, msm\t %.3lf, rel. err\t %.3lf" % (i + 1, B_sim[i], B_msm[i], np.abs(B_sim[i] - B_msm[i])/ B_msm[i] * 100 ))


for i in range(i + 1, len(B_sim)):
    print("measurement %d: sim\t %.3lf" % (i + 1, B_sim[i]))
    
    
# -------------------------------- create a plot ------------------------------
# plot magnetic flux density over position
plt.figure(1)
plt.clf()
xi = range(0, 25)
plt.plot([7, 7], [0, 1.6], "--k")
plt.plot([18, 18], [0, 1.6], "--k")
plt.plot(xi, B_msm[0:25], 'o', label="measured")
plt.plot(xi, B_sim[0:25], '*', label="calculated")
plt.title("measured and calculated values")
plt.ylabel("B")
plt.xlabel("z/x/z - Position")
plt.ylim(0, 1.6)
plt.grid()

plt.annotate("air gap",xy=(7, 1), arrowprops=dict(arrowstyle='->'), xytext=(10, 1.2))
plt.annotate("corner",xy=(18, 0.2), arrowprops=dict(arrowstyle='->'), xytext=(15, 0.6))

plt.legend(loc=0)
plt.show(block=False)
    
loadView(filename="./viewOpt.txt")
cmdInput(locals(), globals())

