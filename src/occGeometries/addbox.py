from netgen.occ import *
import netgen.meshing as meshing
from ngsolve import Mesh

geo = OCCGeometry("Induktivitaet.step")
# geo.Glue()
shape = geo.shape
# shape = Glue(geo.shape)

subsol = shape.SubShapes(TopAbs_ShapeEnum.SOLID)
shape = Glue(subsol)


box = Box( (-20,-25, -5), (20, 20, 25) )

coil = subsol[1]+subsol[2]+subsol[3]

air = box-coil-subsol[0]
# air = box-shape
air.mat("air")



subsol.append (air)

# newgeo = OCCGeometry(subsol)
# newgeo = OCCGeometry([shape,air])
newgeo = OCCGeometry([shape, air])

mesh = Mesh(newgeo.GenerateMesh(maxh=2)) # , perfstepsend=meshing.MeshingStep.MESHEDGES))

Draw (mesh)
print (mesh.GetMaterials())

