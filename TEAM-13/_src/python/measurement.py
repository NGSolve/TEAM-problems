import numpy as np
import matplotlib.pyplot as plt 

import scipy.integrate 

import time

def measure1to25(B, mesh, draw=False):

    Phi_msm = np.zeros(25)
    B_msm = np.zeros(25)

    
    scale = 1e-3
    # +++ measurement 1 - 7 , middle(vertical) sheet +++
    B_norm = B[2].Norm()
    Nx = 100
    Ny = 300
    xi = np.linspace(0, 1.6, Nx) * scale
    yi = np.linspace(-25, 25, Ny) * scale
    zi = np.array([0, 10, 20, 30, 40, 50, 60]) * scale

    A_sheet = 50 * 1.6 * scale**2

    Int_val = np.zeros([len(xi), len(yi)])



    for k in range(len(zi)):
        print("\rmeasurement: " + str(k + 1), end="")
        
        # measure flux
        for i in range(len(xi)):
            for j in range(len(yi)):
                    Int_val[i, j] = float(B_norm(mesh(xi[i], yi[j], zi[k])))
       

         
        Phi_msm[k] =  scipy.integrate.simps(scipy.integrate.simps(Int_val, yi), xi) 
        B_msm[k] = Phi_msm[k] / A_sheet
        
        if draw:
            plt.figure(1)
            X, Y = np.meshgrid(xi, yi)
            plt.clf()
            plt.contourf(X, Y, np.transpose(Int_val), vmin=0, vmax=1.6)
            plt.colorbar()
            plt.axis('equal')
            plt.title(str(k+1))
            plt.draw()
            plt.show(block=False)   
            
            time.sleep(0.5)
        
     # +++ measurement 8 - 18 , back right (vertical) shpeet +++
    B_norm = B[0].Norm()
    Nx = 300
    Ny = 100
    xi = np.array([2.1, 10, 20, 30, 40, 50, 60, 80, 100, 110, 122.1]) * scale
    yi = np.linspace(15, 65, Nx) * scale
    zi = np.linspace(60, 63.2, Ny) * scale

    A_sheet = 50 * 3.2 * scale**2

    Int_val = np.zeros([len(yi), len(zi)])

    for k in range(len(xi)):
        print("\rmeasurement: " + str(k + 7 + 1), end="")
        # measure flux
        for i in range(len(yi)):
            for j in range(len(zi)):
                    Int_val[i, j] = float(B_norm(mesh(xi[k], yi[i], zi[j])))
       

         
        Phi_msm[k + 7] =  scipy.integrate.simps(scipy.integrate.simps(Int_val, zi), yi) 
        B_msm[k + 7] = Phi_msm[k + 7] / A_sheet
        if draw:
            X, Y = np.meshgrid(yi, zi)
            plt.clf()
            plt.contourf(X, Y, np.transpose(Int_val), vmin=0, vmax=1.6)
            plt.colorbar()
            plt.axis('equal')
            plt.title(str(k+7+1))
            plt.draw()
            plt.show(block=False)   
            
            time.sleep(0.5)
    # +++ measurement 8 - 18 , back right (vertical) sheet +++
    B_norm = B[2].Norm()    
    Nx = 100
    Ny = 300
    xi = np.linspace(122.1, 125.3, Nx) * scale
    yi = np.linspace(15, 65, Ny) * scale
    zi = np.array([60, 50, 40, 30, 20, 10, 0]) * scale

    A_sheet = 50 * 3.2 * scale**2

    Int_val = np.zeros([len(xi), len(yi)])

    for k in range(len(zi)):
        print("\rmeasurement: " + str(k + 18 + 1), end="")
        # measure flux
        for i in range(len(xi)):
            for j in range(len(yi)):
                    Int_val[i, j] = float(B_norm(mesh(xi[i], yi[j], zi[k])))
       

         
        Phi_msm[k + 18] =  scipy.integrate.simps(scipy.integrate.simps(Int_val, yi), xi) 
        B_msm[k + 18] = Phi_msm[k + 18] / A_sheet

        if draw:
            X, Y = np.meshgrid(xi, yi)
            plt.clf()
            plt.contourf(X, Y, np.transpose(Int_val), vmin=0, vmax=1.6)
            plt.colorbar()
            plt.axis('equal')
            plt.title(str(k+18+1))
            plt.draw()
            plt.show(block=False)   
            
            time.sleep(0.5)    
    
    print("")
    return B_msm #, Phi_msm
    
def measure26to36(B, mesh):
    B_msm = np.zeros(11)
    B_norm = B.Norm()
    scale = 1e-3
    # +++ measurement 1 - 7 , middle(vertical) sheet +++
    xi = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]) * scale
    y = 20*scale
    z = 55*scale


    for k in range(len(xi)):
        print("\rmeasurement: " + str(k + 1 + 25), end="")
        # measure flux
        B_msm[k] = float(B_norm(mesh(xi[k], y, z)))
    print("")
    return B_msm


def measure37to40(B, mesh):
    B_msm = np.zeros(4)
    B_norm = B.Norm()
    scale = 1e-3
    # +++ measurement 1 - 7 , middle(vertical) sheet +++
    xi = np.array([2.2, 2.0, 1.5, 1.5]) * scale
    yi = np.array([15.1, 14.9, 0.0, 0.0]) * scale
    zi = np.array([60.1, 50.9, 55.0, 25.0]) * scale


    for k in range(len(xi)):
        print("\rmeasurement: " + str(k + 1 + 36), end="")
        # measure flux
        B_msm[k] = float(B_norm(mesh(xi[k], yi[k], zi[k])))
    print("")
    return B_msm
    
    
