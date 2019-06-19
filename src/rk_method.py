import numpy as np


class TimeIntegrator(object):
    def __init__(self, flux_calculator):
        self.flux_calculator = flux_calculator


class RK3Integrator(TimeIntegrator):
    def timestep(G):

        #Compute initial face fluxes
        face_fluxes = self.flux_calculator.compute_face_fluxes(G) 
        
        #RK step # 1
        G1 = []
        k1 = []
        for i, cell in enumerate(G):
            #DEBUG
            if Debug_Flag == 1:
                print 'Computing Fluxes for Cell: ' + str(i+1)
            
            #First RK Step Variable
            for var in cell.keys():
                k1.append( (face_fluxes[i]['left'] - face_fluxes[i]['right']) * (1.0 / dx) )
            
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










