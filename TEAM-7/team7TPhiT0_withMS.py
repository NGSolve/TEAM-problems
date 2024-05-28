from netgen.occ import *
from ngsolve import *

import matplotlib.pyplot as plt
import numpy as np
from netgen.meshing import MeshingParameters
import netgen.gui


from myPackage import assert_almost

xi_msm = np.array([0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]) * 1e-3

def generateMesh(diff = 0.019/2):


    plate = Box(Pnt(0, 0, -diff), Pnt(0.294, 0.294, 0.019-diff)) 
    plate.name="plate"
    plate.faces.name="plate_coat"
    plate.faces.Max(Z).name="plate_top"
    plate.faces.Min(Z).name="plate_bottom"

    hole = Box(Pnt(0.018, 0.018, -diff), Pnt(0.126, 0.126, 0.019 - diff))
    hole.name="hole"
    hole.faces.name="interface"
    hole.faces.Max(Z).name="hole"
    hole.faces.Min(Z).name="hole"


    plate -= hole
    plate.faces.maxh=0.01
    plate.faces.col=(0.3,0.3,0.3)



    wp = WorkPlane(Axes(p=(0.294,0.000,0.049-diff), n=Z, h=Y)) # origin in 
    coil2DCorner = (wp.Rectangle(0.050, 0.050).Face() * wp.MoveTo(0.050, 0.050).Circle(0.050).Face() - wp.Circle(0.025).Face()) 
    coil2DLimb = wp.MoveTo(0.050, 0).Rectangle(0.100, 0.025).Face()
    coil2DCorner.faces.name = "coil_corner"
    coil2DLimb.faces.name = "coil_limb"



    coil2Dpart = Glue([coil2DCorner, coil2DLimb])
    ax = Axis(p=(0.194,0.100,0.049-diff), d=(0, 0, 1))
    coil2D = Glue([coil2Dpart.Rotate(ax, ang) for ang in [0, 90, 180, 270]])

    coil2DInner = wp.MoveTo(0.025, 0.025).Rectangle(0.150, 0.150).Face()
    coil2DInner -= coil2D
    coil2DInner.faces.name ="coil_inner"

    coil2D = Glue([coil2D, coil2DInner])

    coil = coil2D.Extrude(0.100)

    for i in range(len(coil.faces)):
        coil.faces[i].name = "coil"
    #Draw(coil)

    bounding_box = Box(Pnt(-.2, -.2, -.2-diff), Pnt(0.5, 0.5, 0.5-diff))
    bounding_box.faces.Max(X).name="right"
    bounding_box.faces.Min(X).name="left"
    bounding_box.faces.Max(Y).name="back"
    bounding_box.faces.Min(Y).name="front"
    bounding_box.faces.Max(Z).name="top"
    bounding_box.faces.Min(Z).name="bottom"
    bounding_box.name="air"

    bounding_box = bounding_box-coil- plate-hole


    geo = OCCGeometry(Glue([bounding_box, coil, plate, hole]))
    # geo = OCCGeometry(Glue([coil]))

    mp = MeshingParameters ( maxh =1)
    if True:
        # for A1B1 measurement
        [mp.RestrictH(x, 0.072, 0.034-diff, h= 0.01/3) for x in xi_msm];
        # for A2B2 measurement
        [mp.RestrictH(x, 0.144, 0.034-diff, h= 0.01/3) for x in xi_msm];
        # for A3B3 measurement
        [mp.RestrictH(x, 0.072, 0.019-diff, h= 0.01/3) for x in xi_msm];
        # for A4B4 measurement
        [mp.RestrictH(x, 0.072, 0.000-diff, h= 0.01/3) for x in xi_msm];

    mesh = Mesh(geo.GenerateMesh(maxh=0.1, mp=mp))
    
    print("mats", set(mesh.GetMaterials()))
    print("bnd", set(mesh.GetBoundaries()))
    print("number of elements", mesh.ne)

    return mesh


def getT0(mesh):
    val_inner = 1
    val_outer = 0

    center_x, center_y = 0.294-0.100, 0.100
    d = 0.025


    CF_limb = CF(0)
    #right limb
    xs = x - center_x - 0.075
    CF_limb = IfPos(xs +0.01,  -(xs - d) / d, 0)
    # left limb
    xs = x - center_x + 0.075
    CF_limb += IfPos(xs -0.01, 0, (xs + d)/d)
    # top limb
    ys = y - center_y - 0.075
    CF_limb += IfPos(ys +0.01,  -(ys - d) / d, 0)
    # bottom limb
    ys = y - center_y + 0.075
    CF_limb += IfPos(ys -0.01, 0, (ys + d)/d)


    CF_corner = CF(0)

    # bottom right
    corner_x, corner_y  = center_x + 0.05, center_y-0.05
    r = sqrt((x - corner_x)**2 + (y - corner_y)**2)
    CF_corner += IfPos(x - center_x, IfPos(y - center_y, 0, - (r - 0.05)/d ), 0)

    #top right
    corner_x, corner_y  = center_x + 0.05, center_y+0.05
    r = sqrt((x - corner_x)**2 + (y - corner_y)**2)
    CF_corner += IfPos(x - center_x, IfPos(y - center_y, - (r - 0.05)/d,0  ), 0)

    # bottom left
    corner_x, corner_y  = center_x - 0.05, center_y-0.05
    r = sqrt((x - corner_x)**2 + (y - corner_y)**2)
    CF_corner += IfPos(x - center_x, 0, IfPos(y - center_y, 0, - (r - 0.05)/d ))

    #top right
    corner_x, corner_y  = center_x - 0.05, center_y+0.05
    r = sqrt((x - corner_x)**2 + (y - corner_y)**2)
    CF_corner += IfPos(x - center_x, 0, IfPos(y - center_y, - (r - 0.05)/d,0  ))
    

    
    T0 = mesh.MaterialCF({"coil_inner": val_inner, "coil_limb":CF_limb, "coil_corner":CF_corner}, default=val_outer)

    return T0 * CF((0, 0, 1))

def calcWithT0(mesh, T0, order0 = 0):
    print("-"*50)
    print("with T0")
    print("-"*50)
    mu = 4e-7*np.pi
    sigma = mesh.MaterialCF({"plate":3.526e7}, default=1)
    rho = 1/sigma

    freq = Parameter(50)
    J0 = 2742 / (0.025 * 0.100)
    print(f"J0 = {J0}  A/m**2")
    omega = freq * np.pi * 2
    
    

    fesT = HCurl(mesh, order=order0, dirichlet="plate.*|hole", complex=True, definedon="plate|hole")
    fesPhi = H1(mesh, order=order0+1, dirichlet="top|back|right|bottom", complex=True)
    fes = FESpace([fesPhi, fesT])
    
    print("free dofs", sum(fes.FreeDofs()))
    trials, tests = fes.TnT()
    uPhi, vPhi = trials[0], tests[0]
    uT, vT = trials[1], tests[1]

    sol = GridFunction(fes, "sol")
    Phi = sol.components[0]
    T = sol.components[1]
    
    J = curl(T)
    H = T - grad(Phi)

    B = mu * H
    E = rho * J
    
    a = BilinearForm(fes, symmetric=True)
    a += rho * curl(uT) * curl(vT) * dx("hole|plate")
    a += 1j * omega * mu * (uT - grad(uPhi)) * (vT - grad(vPhi)) * dx

    a += 1e-1 * uPhi * vPhi * dx("hole|plate")


    
    
    f = LinearForm(fes)
    # f += -1j*omega * mu * J0 * J_coil* curl(vT) * dx("coil.*")
    # f += -1j*omega * mu * Cross(n, J0*J_coil_bnd) * (vT.Trace() - grad(vPhi).Trace())* ds("coil")
    f += -1j * omega  * mu * J0 * 0.025 * T0 * (vT - grad(vPhi)) * dx("coil.*")

    

    f.Assemble()
    print(sum(f.vec))
    
    pre = Preconditioner(a, type="direct", inverse="sparsecholesky")
    print("start")
    with TaskManager():
        solvers.BVP(bf=a, lf = f, pre=pre, gf=sol, tol=1e-20)
        

    return {"sol":sol, "H": H, "B":B, "J":J}

def getMsmValues():
    Bz_A1B1_50Hz_ref = np.array([ -4.9  -1.16j, -17.88 +2.48j, -22.13 +4.15j, -20.19 +4.j,   -15.67 +3.07j,
        0.36 +2.31j,  43.64 +1.89j,  78.11 +4.97j,  71.55+12.61j , 60.44+14.15j,
        53.91+13.04j,  52.62+12.4j,  53.81+12.05j, 56.91+12.27j,  59.24+12.66j,
        52.78 +9.96j,  27.61 +2.26j])
                
    Bz_A2B2_50Hz_ref = np.array([-1.83-1.63j, -8.5-0.6j, -13.6-0.43j, -15.21+0.11j, -14.48+1.26j, -5.62+3.4j,
        28.77+6.53j, 60.34+10.25j, 61.84+11.83j, 56.64+11.83j, 53.4+11.01j, 52.36+10.58j, 53.93+10.8j, 56.82+10.54j, 
        59.48+10.62j, 52.08+9.03j, 26.56+1.79j])

    Jy_A3B3_50Hz_ref = np.array([0.249-0.629j,  0.685-0.873j,  0+0j,  0+0j,  0+0j,  0+0j,  0+0j,  -0.015-0.593j,
        -0.103-0.249j,  -0.061-0.101j,  -0.004-0.001j,  0.051+0.087j,  0.095+0.182j,  0.135+0.322j,  
        0.104+0.555j,  -0.321+0.822j,  -0.687+0.855j,])


    Jy_A4B4_50Hz_ref = np.array([0.461-0.662j,  0.621-0.664j,  0+0j,  0+0j,  0+0j,  0+0j,  0+0j,  1.573-1.027j,
        0.556-0.757j,  0.237-0.364j,  0.097-0.149j,  -0.034+0.015j,  -0.157+0.154j,  -0.305+0.311j,
        -0.478+0.508j,  -0.66+0.747j,  -1.217+1.034j])
    
    return Bz_A1B1_50Hz_ref, Bz_A2B2_50Hz_ref, Jy_A3B3_50Hz_ref, Jy_A4B4_50Hz_ref




def calcWithT0MS(mesh, T0, order0 = 0):
    print("-"*50)
    print("with T0 MS")
    print("-"*50)
    mu = 4e-7*np.pi

    sigmaPlate = 3.526e7
    sigma = mesh.MaterialCF({"plate":sigmaPlate}, default=1)
    rho = 1/sigma

    freq = Parameter(50)
    J0 = 2742 / (0.025 * 0.100)
    print(f"J0 = {J0}  A/m**2")
    omega = freq * np.pi * 2


    from MS_helper_functions import cl_Phi, cl_curlcurlMS, cl_gradgradMS
    cl_Phi.numSheets = 1
    cl_Phi.dFe = 0.019
    cl_Phi.d0 = 0.01
    cl_Phi.mesh = mesh

    cl_Phi.modelHalfAir = False
    cl_Phi.orientation = 2 

    dom = "plate|hole"
    orderPhi = [
        cl_Phi(1, fes_order=1, material=dom, dirichlet="hole|plate_top|plate_bottom"), 
        cl_Phi(2, fes_order=1, material=dom, dirichlet="hole|plate_top|plate_bottom", inAir=False), 
    ]

    orderT = [
        cl_Phi(0, fes_order=1, material=dom, dirichlet="hole|plate.*", inAir=False, modelHalfAir=False, nograds=True, useAbsolutes=False),
    ]


    VSpace = []
    VSpace.append(H1(mesh, order=order0+1, dirichlet="top|back|right|bottom", complex=True)) 
        
    # ui * phi i
    for phi_i in orderPhi: 
        VSpace.append(H1(mesh, order=phi_i.fes_oder+1, definedon=phi_i.material, dirichlet=phi_i.dirichlet, complex=True))

    for phi_i in orderT: 
        VSpace.append(HCurl(mesh, order=phi_i.fes_oder, definedon=phi_i.material, dirichlet=phi_i.dirichlet, complex=True, nograds=phi_i.nograds))

    VSpace = FESpace(VSpace)


    # multiscale container
    sol = GridFunction(VSpace, "sol")
    print("MS ndofs", sum(VSpace.FreeDofs()))
    gradgradMS = cl_gradgradMS(orderPhi, sol, addPhi0Outer=True, secondOrder=False)
    curlcurlMS = cl_curlcurlMS(orderT, sol, eddy_inplane=False, istart = len(gradgradMS.orderPhi))
    gradgradMS.addCurlCurlMS(curlcurlMS)
    curlcurlMS.addGradGradMS(gradgradMS)


    # BLF
    a = BilinearForm(VSpace, symmetric=True)
    a += 1j * omega * mu  * grad(gradgradMS.trials[0]) * grad(gradgradMS.tests[0]) * dx("air|coil.*")
    a += 1j * 1e2  * gradgradMS.trials[0] * gradgradMS.tests[0] * dx(dom)
    a += 1j * omega * gradgradMS.getIntegrand4BFI(gradgradMS.gradu_pack, gradgradMS.gradv_pack, mu, mu) * dx(dom)

    a += curlcurlMS.getIntegrand4BFI(curlcurlMS.curlu_pack, curlcurlMS.curlv_pack, 1/sigmaPlate, 1/1, conductivity=False) * dx("plate")
    a += curlcurlMS.getIntegrand4BFI(curlcurlMS.curlu_pack, curlcurlMS.curlv_pack, 1/1, 1/1, conductivity=False) * dx("hole")
    

    f = LinearForm(VSpace)
    f += -1j * omega * gradgradMS.getIntegrand4LFI(J0 * 0.025 * T0, gradgradMS.gradv_pack[:], mu, mu) * dx("coil.*")
        
    print("solving")
    with TaskManager():    
        prec = Preconditioner(a,type="direct")  
        solvers.BVP(bf = a, lf= f, pre=prec, gf=sol, maxsteps=50, tol = 1e-25)

    H = sum(gradgradMS.gradsol_comp)
    J = sum(curlcurlMS.curlsol_comp)
                
    return {"sol":sol, "H": H, "B":mu*H, "J":J}
def run():
    from myPackage import myBreak, drawBndAll, drawDomainsAll
    
    diff = 0.019/2
    try:
        mesh = Mesh("mesh.vol")
    except:
        
        mesh = generateMesh(diff=diff)
        mesh.ngmesh.Save("mesh.vol")

    mesh.Curve(3)

    Vol_coil_should = 0.100 * (0.025*0.100*4 + (0.050**2-0.025**2) * np.pi)
    Vol_plate_should = 0.019 * (0.294*0.294-0.108*0.108)
    assert_almost(Integrate(1, mesh, definedon=mesh.Materials("coil_limb|coil_corner")),  Vol_coil_should, eps=1e-4)
    assert_almost(Integrate(1, mesh, definedon=mesh.Materials("plate.*")),  Vol_plate_should, eps=1e-4)
        
    T0 = getT0(mesh)


    retMS = calcWithT0MS(mesh, T0)
    BMS = retMS["B"]
    JMS = retMS["J"]

    Draw(BMS, mesh, "BMS")
    Draw(JMS, mesh, "JMS")


    ret = calcWithT0(mesh, T0)

    B = ret["B"]
    J = ret["J"]

    Draw(B, mesh, "B")
    Draw(J, mesh, "J")


    Bz_A1B1_50Hz_ref, Bz_A2B2_50Hz_ref, Jy_A3B3_50Hz_ref, Jy_A4B4_50Hz_ref = getMsmValues()


    Bz_A1B1 = np.array([(B[2])(mesh(x, 0.072, 0.034- diff)) for x in xi_msm])
    Bz_A1B1MS = np.array([(BMS[2])(mesh(x, 0.072, 0.034 - diff)) for x in xi_msm])
    print(Bz_A1B1)
    plt.figure(1)
    plt.clf()
    plt.title("Evaluate B.z on Line A1 - B1")
    plt.plot(xi_msm, 1e4 * Bz_A1B1.real, "b-x", label="sim Bz.real")
    plt.plot(xi_msm, -1e4 * Bz_A1B1.imag, "r-x", label="sim Bz.imag")
    plt.plot(xi_msm, 1e4 * Bz_A1B1MS.real, "b-o", label="sim Bz.real")
    plt.plot(xi_msm, -1e4 * Bz_A1B1MS.imag, "r-o", label="sim Bz.imag")
    plt.plot(xi_msm, Bz_A1B1_50Hz_ref.real, "b--x", label="ref Bz.real")
    plt.plot(xi_msm, Bz_A1B1_50Hz_ref.imag, "r--x", label="ref Bz.imag")
    plt.ylabel("B in Gauss = $10^{-4}$ T")
    plt.legend()
    plt.show()


    Jy_A3B3 = np.array([J[1](mesh(x, 0.072, 0.019-1e-5)) for x in xi_msm])
    # Jy_A3B3MS = np.array([JMS[1](mesh(x, 0.072, 0.019-1e-5)) for x in xi_msm])

    plt.figure(3)
    plt.clf()
    plt.title("Evaluate J.y on Line A3 - B3")
    plt.plot( xi_msm, 1e-6 *Jy_A3B3.real, "b-x", label="sim Jy.real")
    plt.plot( xi_msm, 1e-6 * Jy_A3B3.imag, "r-x", label="sim Jy.imag")
    # plt.plot( xi_msm, 1e-6 *Jy_A3B3MS.real, "b-o", label="sim Jy.real")
    # plt.plot( xi_msm, 1e-6 * Jy_A3B3MS.imag, "r-o", label="sim Jy.imag")
    plt.plot(xi_msm, Jy_A3B3_50Hz_ref.real, "b--x", label="ref Jy.real")
    plt.plot(xi_msm, Jy_A3B3_50Hz_ref.imag, "r--x", label="ref Jy.imag")
    plt.ylabel("J in $10^{6}$ A/mm$^2$")
    plt.ylim([min(list(Jy_A3B3_50Hz_ref.real)+list(Jy_A3B3_50Hz_ref.imag)), max(list(Jy_A3B3_50Hz_ref.real)+list(Jy_A3B3_50Hz_ref.imag))])
    plt.legend()
    plt.show()


    myBreak(locals(), globals(), file=__file__.split('/')[-1])