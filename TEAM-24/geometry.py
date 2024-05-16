from netgen.csg import *
from netgen.meshing import *

def MakeGeometry():

    topbot = Plane (Pnt(0, 0,-0.0127), Vec(0,0,-1)) * Plane (Pnt(0,0,0.0127), Vec(0,0,1))

    stator = Cylinder (Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.1045) - Cylinder(Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.0831)
    stator = stator + OrthoBrick (Pnt(-0.0139,-0.090,-0.020), Pnt(0.0139, 0.090, 0.020)) - Cylinder(Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.05375).bc("integrate")
    stator = stator * topbot

    SetTransformation(dir=3, angle=22)
    p1 = Pnt(-0.020, -0.080, -0.0127)
    p2 = Pnt(0.020, 0.080, 0.0127)
    rotor = Cylinder (Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.0341)  + Plane(p1, Vec(-1,0,0))*Plane(p1, Vec(0,-1,0))*Plane(p1, Vec(0,0,-1))*Plane(p2, Vec(1,0,0))*Plane(p2, Vec(0,1,0))*Plane(p2, Vec(0,0,1))
    rotor = rotor * Cylinder(Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.05105).bc("rotor_end") - Cylinder(Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.0254)
    rotor = rotor * topbot
    SetTransformation()

    c1 = OrthoBrick(Pnt(-0.041, -0.100, -0.0155), Pnt(-0.017, 0.100, 0.0155))
    c2 = Cylinder (Pnt(-0.017,-0.100, 0.0155), Pnt(-0.017,0.100,0.0155), 0.024) * Plane(Pnt(-0.017,-0.100,0.0155), Vec(1,0,0)) *  Plane(Pnt(-0.017,-0.100,0.0155), Vec(0,0,-1))
    c3 = OrthoBrick(Pnt(-0.017, -0.100, 0.0155), Pnt(0.017, 0.100, 0.0395))
    c4 = Cylinder (Pnt(0.017,-0.100,0.0155), Pnt(0.017,0.100,0.0155), 0.024) * Plane(Pnt(0.017,-0.100,0.0155), Vec(-1,0,0)) *  Plane(Pnt(0.017,-0.100,0.0155), Vec(0,0,-1))    
    c5 = OrthoBrick(Pnt(0.017, -0.100, -0.0155), Pnt(0.041, 0.100, 0.0155))
    c6 = Cylinder (Pnt(0.017,-0.100,-0.0155), Pnt(0.017,0.100,-0.0155), 0.024) * Plane(Pnt(0.017,-0.100,-0.0155), Vec(-1,0,0)) *  Plane(Pnt(0.017,-0.100,-0.0155), Vec(0,0,1))
    c7 = OrthoBrick(Pnt(-0.017, -0.100, -0.0395), Pnt(0.017, 0.100, -0.0155))
    c8 = Cylinder (Pnt(-0.017,-0.100,-0.0155), Pnt(-0.017,0.100,-0.0155), 0.024) * Plane(Pnt(-0.017,-0.100,-0.0155), Vec(1,0,0)) *  Plane(Pnt(-0.017,-0.100,-0.0155), Vec(0,0,1))        
    coilcatA = OrthoBrick(Pnt(-0.100,-0.0715, -0.080), Pnt(0.100,-0.0545,0.080)).maxh(0.010)
    coilcatB = OrthoBrick(Pnt(-0.100,0.0545, -0.080), Pnt(0.100,0.0715,0.080)).maxh(0.010)

    box =  Cylinder (Pnt(0,0,-0.020), Pnt(0,0,0.020), 0.1505) * Plane(Pnt(0,0,-0.100),Vec(0,0,-1)) * Plane(Pnt(0,0,0.100),Vec(0,0,1))
    
    sym = Plane (Pnt(0,0,0), Vec(0,0,1)).bc("sym")
    geo = CSGeometry()

    objmaxh = 0.005
    geo.Add ((sym*stator.maxh(objmaxh)).mat("stator"), maxh=objmaxh)
    geo.Add ((sym*rotor.maxh(objmaxh)).mat("rotor"), maxh=objmaxh)
    geo.Add ((sym*c1*coilcatA).mat("coil1a"), maxh=objmaxh)
    geo.Add ((sym*c2*coilcatA).mat("coil2a"), maxh=objmaxh)
    geo.Add ((sym*c3*coilcatA).mat("coil3a"), maxh=objmaxh)
    geo.Add ((sym*c4*coilcatA).mat("coil4a"), maxh=objmaxh)
    geo.Add ((sym*c5*coilcatA).mat("coil5a"), maxh=objmaxh)
    geo.Add ((sym*c6*coilcatA).mat("coil6a"), maxh=objmaxh)
    geo.Add ((sym*c7*coilcatA).mat("coil7a"), maxh=objmaxh)
    geo.Add ((sym*c8*coilcatA).mat("coil8a"), maxh=objmaxh)
    geo.Add ((sym*c1*coilcatB).mat("coil1b"), maxh=objmaxh)
    geo.Add ((sym*c2*coilcatB).mat("coil2b"), maxh=objmaxh)
    geo.Add ((sym*c3*coilcatB).mat("coil3b"), maxh=objmaxh)
    geo.Add ((sym*c4*coilcatB).mat("coil4b"), maxh=objmaxh)
    geo.Add ((sym*c5*coilcatB).mat("coil5b"), maxh=objmaxh)
    geo.Add ((sym*c6*coilcatB).mat("coil6b"), maxh=objmaxh)
    geo.Add ((sym*c7*coilcatB).mat("coil7b"), maxh=objmaxh)
    geo.Add ((sym*c8*coilcatB).mat("coil8b"), maxh=objmaxh)
    geo.Add ((sym*(box-stator-rotor- (c1+c2+c3+c4+c5+c6+c7+c8)*(coilcatA+coilcatB))).mat("air"))
    
    return geo


if __name__ == "__main__":
    geo = MakeGeometry()
    import ngsolve
    mesh = ngsolve.Mesh(geo.GenerateMesh())
    ngsolve.Draw(mesh)
