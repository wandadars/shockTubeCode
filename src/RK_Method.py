#######################################################################################
# The purpose of this function is to use a Runge-kutta time stepping scheme to advance
# the conserved variables in time.
#
# Author: Christopher Neal
#
# Date: 06-18-2015
# Updated: 06-24-2015
#
#######################################################################################
#
# Information: This function takes in the conserved variables for all cells and computes
#          the new value of the conserved variable at the new time using RK method.


import numpy as np
import Compute_Face_Fluxes as FF
from Input_Params import *

class TimeIntegrator(object):
    def __init__(self):
        pass

class RK3Integrator(TimeIntegrator):
    def timestep(G):
        #Compute initial face fluxes
        FaceFluxes = FF.compute_cell_face_fluxes(G)
        
        #RK step # 1
        G1 = np.zeros((Nx,3))
        k1 = np.zeros((Nx,3))

        for i in range(0,Nx,1):
            #DEBUG
            if(Debug_Flag == 1):
                print 'Computing Fluxes for Cell: ' + str(i+1)
            
            #First RK Step Variable
            k1[i,:] = (FaceFluxes[i,:,0] - FaceFluxes[i,:,1])*(1.0/dx)
            G1[i,:] = G[i,:] + dt*(k1[i,:]) # Explicit Euler approximation


        #Compute fluxes for Second RK Variable
        FaceFluxes = FF.compute_cell_face_fluxes(G1)

        #RK step # 2
        G2 = np.zeros((Nx,3))
        k2 = np.zeros((Nx,3))

        for i in range(0,Nx,1):
            #Second RK Step Variable
            k2[i,:] = (FaceFluxes[i,:,0] - FaceFluxes[i,:,1])*(1.0/dx)
            G2[i,:] = G[i,:] + (dt/2.0)*(0.5)*(k1[i,:] + k2[i,:])

        #Compute fluxes for Third RK Variable
        FaceFluxes = FF.compute_cell_face_fluxes(G2)

        #RK Step # 3
        k3 = np.zeros((Nx,3))
        for i in range(0,Nx,1):
            #Third RK Step Variable
            k3[i,:] = (FaceFluxes[i,:,0] - FaceFluxes[i,:,1])*(1.0/dx)


        G_Final = G + (1.0/6.0)*dt*(k1 + k2 + 4.0*k3)

        return G_Final










