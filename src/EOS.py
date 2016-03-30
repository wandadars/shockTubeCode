#######################################################################################
# The purpose of this function is to compute a thermodynamic property from an IDEAL 
# equation of state and return the result.
#
# Author: Christopher Neal
#
# Date: 06-18-2015
# Updated: 06-18-2015
#
#######################################################################################
#
# Information: This function takes in the conserved variables for a cell and computes

def IDEAL_EOS_P(CellValue):
	
	from Input_Params import *

	# CellValues is a 3 element vector that contains the conserved variables that are
	# used for computing a value from the equation of state.

	#This particular function computes pressure
	# P = (gamma-1)*(rho*e - (1/2)*rho*u*u)
	P = (gamma-1.0)*(CellValue[2] - 0.5*CellValue[1]*(CellValue[1]/CellValue[0]) )

	return P

