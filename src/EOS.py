"""
The purpose of this function is to compute a thermodynamic property from an IDEAL 
equation of state and return the result.
  Information: This function takes in the conserved variables for a cell and computes
"""

class EOS(object):
    def __init__(self):
        pass


class IdealEOS(EOS):
    def __init__(self, gamma, r_gas):
        self.gamma = gamma
        self.r_gas = r_gas

    def energy(pressure, density):
        return pressure / (density * (self.gamma - 1))

    def pressure(density, energy, velocity):
        return (self.gamma - 1.0) * (density*energy - 0.5*velocity**2))

    def density(pressure, temperature)
        return pressure / (self.r_gas * temperature)

def IDEAL_EOS_P(CellValue):
	from Input_Params import *

	# CellValues is a 3 element vector that contains the conserved variables that are
	# used for computing a value from the equation of state.

	#This particular function computes pressure
	# P = (gamma-1)*(rho*e - (1/2)*rho*u*u)
	P = (gamma-1.0)*(CellValue[2] - 0.5*CellValue[1]*(CellValue[1]/CellValue[0]) )

	return P

