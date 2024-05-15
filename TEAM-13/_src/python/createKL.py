import numpy as np

import scipy as sp
from scipy import interpolate



def __getH(B):
    B_KL = np.array([0, 25e-4, 50e-4, 125e-4, 250e-4, 5e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8])
    H_KL = np.array([0, 16, 30, 54, 93, 143, 191, 210, 222, 233, 247, 258, 272, 289, 313, 342, 377, 433, 509, 648, 933, 1228, 1934, 2913, 4993, 7189, 9423])

    mu0 = 4e-7*np.pi

    a = -2.381e-10
    b = 2.327e-5
    c = 1.590

    Ms = 2.16

    B = np.array(B)
    H = np.zeros(np.shape(B))
    ind_interpolate = (B <= 1.8)
    ind_saturation = (B > 2.22)
    ind_quad = np.logical_not(np.logical_or(ind_saturation, ind_interpolate))

    H[ind_quad] = ((mu0**2+2*b*mu0-4*a*c+b**2+4*B[ind_quad]*a)**(1/2)-mu0-b)/(2*a)
    H[ind_interpolate] = sp.interpolate.interp1d(B_KL, H_KL, kind='cubic')(B[ind_interpolate])
    H[ind_saturation] = (B[ind_saturation] - Ms) / mu0

    return H


def createKL():

    B_KL = np.array([0, 25e-4, 50e-4, 125e-4, 250e-4, 5e-2, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8])
    H_KL = np.array([0, 16, 30, 54, 93, 143, 191, 210, 222, 233, 247, 258, 272, 289, 313, 342, 377, 433, 509, 648, 933, 1228, 1934, 2913, 4993, 7189, 9423])
    Bi = np.hstack([B_KL, np.linspace(max(B_KL), 5, 50)])

    Hi = __getH(Bi)

    return Hi, Bi
    

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    H_KL, B_KL = createKL()
    
    plt.figure(1)
    plt.clf()
    plt.plot(H_KL, B_KL, '.-r')
    plt.xlim(0, 2000)
    plt.ylim(0, 2.12)
    plt.grid()
    plt.title("B/H-curve")
    plt.xlabel("H / A/m")
    plt.ylabel("B / T")
    plt.show(block=False)
    
    print("H of curve:")
    print(H_KL)

    print("B of curve:")
    print(B_KL)
    
    try:
        from myPackage import cmdInput
        cmdInput(locals(), globals())
    except:
        input("press enter to finish")




