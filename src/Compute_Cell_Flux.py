
import numpy as np
from Input_Params import *
import EOS as eos

def Compute_Cell_Fluxes(CellValue):
    """
    The purpose of this function is to compute the value of the flux variables at the cell
    centers using the conserved variables for the cell.

    Information: This function takes in a cell data 1x3 array that holds the 3 conserved
    variables.
        Output: A 1x3 array of the fluxes for mass, momentum,and energy
    """

	#Allocate the array to hold the flux data
	Cell_Flux = np.zeros(3)

	#Store flux for conservation of mass equation i.e. rho*u
	Cell_Flux[0] = CellValue[1]

	#Store flux for conservation of momentum equation i.e. rho*u^2 + p
	Cell_Flux[1] = (CellValue[1]**2)/CellValue[0] + eos.IDEAL_EOS_P(CellValue)

	#Store flux for conservation of energy equation i.e. u*(rho*e + p)
	Cell_Flux[2] = (CellValue[1]/CellValue[0])*(CellValue[2] + eos.IDEAL_EOS_P(CellValue) )

	return Cell_Flux
