from netgen.csg import *
from ngsolve import *
import numpy as np

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
print ("have energy")


import matplotlib.pyplot as plt

bpts = np.linspace(0,2.2,100)
hpts = [bhcurve(b) for b in bpts]
epts = [energy(b) for b in bpts]

print ("bpts = ", bpts)
print ("hpts = ", hpts)
print ("epts = ", epts)

if True:
    plt.axis([0,200,0,3])
    plt.xlabel("H")
    plt.ylabel("B")
    plt.plot (hpts,bpts)
else:
    plt.xlabel("B")
    plt.ylabel("E")
    plt.plot (bpts, epts)
plt.show()
