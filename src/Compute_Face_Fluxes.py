#######################################################################################
# The purpose of this function is to compute the values of the cell face fluxes the
# hllc method.
#
# Author: Christopher Neal
#
# Date: 06-18-2015
# Updated: 06-18-2015
#
#######################################################################################
#
# Information: This function takes in a cell data Nx X 3 array that holds the 3 conserved
# variables for all cells in the domain.
#
#       Output: A Nx X 3 X 2 array of the fluxes for mass, momentum,and energy at the west
#		and east faces for each cell in the domain

def compute_cell_face_fluxes(G):

	from Input_Params import *
	import numpy as np
        import Face_Flux_Schemes as FFS
	import Compute_BC_Conserved as BCONSERV

	#Allocate array to hold the conserved variables for the left and right cells that surround
	#the face.
	temp = np.zeros((2,3))


        #Compute fluxes at faces by looping over all cells and using a flux function

	#Create array to hold face flux values
	# the 2 is for the east and west fluxes for each cell.
	FaceFluxes = np.zeros((Nx,3,2)) # the 3 is because there are 3 conserved variables


        #The First Cell has a prescribed flux on its left face from the boundary
	#DEBUG
	if(Debug_Flag == 1):
		print 'Debug Info for Cell: '+ str(1) +'\n' 

        #Compute and store left face flux for the first cell
	temp[0,:] = BCONSERV.compute_BC_conserved(0)
	temp[1,:] = G[0,:]

	if(FluxSchemeInput == 0):       #HLLC Flux
                FaceFluxes[0,:,0] = FFS.hllc_flux(temp)

        elif(FluxSchemeInput == 1):     #ROE Flux
                FaceFluxes[0,:,0] = FFS.roe_flux(temp)


        #Compute and store right face flux for the first cell
	 temp[0,:] = G[0,:]
         temp[1,:] = G[1,:]

	if(FluxSchemeInput == 0):	#HLLC Flux
        	FaceFluxes[0,:,1] = FFS.hllc_flux(temp)

	elif(FluxSchemeInput == 1):	#ROE Flux
                FaceFluxes[0,:,1] = FFS.roe_flux(temp)
	


	
        #Inner cells have 2 face fluxes that must be determined
        for j in range(1,Nx-1,1):
		
		#DEBUG
		if(Debug_Flag == 1):
			print 'Debug Info for Cell: '+ str(j+1) +'\n'


		#DEBUG
		if(Debug_Flag == 1):
			print 'Left Cell Face Flux:\n'

		#Compute and store left face flux
		if(FluxSchemeInput == 0):	#HLLC Flux
			#For this type of flux scheme the result of this will be the same as for the previous
			#cell's right face flux
                	FaceFluxes[j,:,0] = FaceFluxes[j-1,:,1]
		elif(FluxSchemeInput == 1):	#Roe Flux
			#For this type of flux scheme the result of this will be the same as for the previous
                        #cell's right face flux
                        FaceFluxes[j,:,0] = FaceFluxes[j-1,:,1]


		#DEBUG
		if(Debug_Flag == 1):
			print 'Right Cell Face Flux:\n'

		#Compute and store right face flux
		temp[0,:] = G[j,:]
                temp[1,:] = G[j+1,:]

		if(FluxSchemeInput == 0):	#HLLC Flux
                	FaceFluxes[j,:,1] = FFS.hllc_flux(temp)
		elif(FluxSchemeInput == 1):	#Roe Flux
                        FaceFluxes[j,:,1] = FFS.roe_flux(temp)




        #The Last Cell has a prescribed flux on its right face from the boundary
	if(Debug_Flag == 1):
		print 'Debug Info for Cell: '+ str(Nx) +'\n'

        
	#DEBUG
	if(Debug_Flag == 1):
	     	print 'Flux across Left Cell Face'

	#Compute and store left face flux
        temp[0,:] = G[Nx-2,:]   #Remember arrays go from 0 to Nx-1 for Nx elements
        temp[1,:] = G[Nx-1,:]

	if(FluxSchemeInput == 0):	#HLLC Flux
        	FaceFluxes[Nx-1,:,0] = FFS.hllc_flux(temp)

	elif(FluxSchemeInput == 1):	#Roe Flux
                FaceFluxes[Nx-1,:,0] = FFS.roe_flux(temp)


	#Compute and store right face flux
        temp[0,:] = G[Nx-1,:]
	temp[1,:] = BCONSERV.compute_BC_conserved(1)

	if(FluxSchemeInput == 0):       #HLLC Flux
                FaceFluxes[Nx-1,:,1] = FFS.hllc_flux(temp)

        elif(FluxSchemeInput == 1):     #Roe Flux
                FaceFluxes[Nx-1,:,1] = FFS.roe_flux(temp)


	return FaceFluxes

