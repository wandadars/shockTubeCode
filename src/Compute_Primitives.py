#######################################################################################
# The purpose of this function is to compute the values of the primitive variables from
# the flow field solution variables.
#
# Author: Christopher Neal
#
# Date: 06-22-2015
# Updated: 06-24-2015
#
#######################################################################################

def compute_primitives(CellValues):

	from Input_Params import *
	import numpy as np

	Primitives = np.zeros((Nx,4)) # rho, u, P, T	

	#Compute Density
	Primitives[:,0] = CellValues[:,0]

	#Compute velocity
	Primitives[:,1] = CellValues[:,1] / CellValues[:,0]

	#Compute Pressure
	Primitives[:,2] = (gamma-1)*(CellValues[:,2] - 0.5*(CellValues[:,1]**2)/CellValues[:,0])
	
	#Compute Temperature: T = ( e - (1/2)u*u) * (gamma-1)/R
	Primitives[:,3] = ( ( CellValues[:,2]/CellValues[:,0] ) - 0.5*(CellValues[:,1]/CellValues[:,0])**2)*(gamma-1)/R_gas


	return Primitives

