#######################################################################################
# This module contains functions for computing important cell state values for different
# flux schems.
#
# Author: Christopher Neal
#
# Date: 09-30-2015
# Updated: 09-30-2015
#
#######################################################################################

from math import sqrt
import numpy as np
import EOS as eos
from Input_Params import *

class CellStates(object):
	"""
	The purpose of this class is to compute a set of state variables for a cell to be
	used in the HLLC flux calculation.
	 Information: This function takes in a 2x3 array of cell left & right cell conservative
	      variables and outputs a 2x5 array that contains the necessary cell variables
	      for the HLLC flux scheme to estimate a flux. A combination of primitive and
	      conserved variables.
	"""
    def __init__(self, adjacent_cell_values):
        self.states = None

    def print_cell_values(self):
        print 'Left Density = ' + str(CellValues[0,0])
        print 'Left Momentum = ' + str(CellValues[0,1])
        print 'Left Energy = ' + str(CellValues[0,2])

        print 'Right Density = ' + str(CellValues[1,0])
        print 'Right Momentum = ' + str(CellValues[1,1])
        print 'Right Energy = ' + str(CellValues[1,2])


class HLLCStates(CellStates):
    def __init__(self):
        self.states = {} 


    def compute_states(self, cell_values):
	#CellValues is a 2x3 array containing the conserved variable information for the
	#left and right states CellValues(0,:)=>Left State & CellValues(1,:)=>Right State

	#States is a 2x5 array that contains the relevant state information for the left
	#and right states that is computed from the conserved cell values i.e. primitive values.
	#the first row is states from the left cell & the second row is right cell states

	#Store the velocities of the states
	states['u_left'] = cell_values[0,1] / cell_values[0,0] #Left velocity
	states['u_right']  cell_values[1,1] / cell_values[1,0] #Right velocity

	#Compute speed of sound(SS)
	states['soundspeed_left'] = np.sqrt( gamma*eos.IDEAL_EOS_P(CellValues[0,:])/CellValues[0,0] ) #Left SS
	states[1,1] = np.sqrt( gamma*eos.IDEAL_EOS_P(CellValues[1,:])/CellValues[1,0] ) #Right SS
	
	#Store density
	states[0,2] = CellValues[0,0]
	states[1,2] = CellValues[1,0]

	#Store Pressure
	states[0,3] = eos.IDEAL_EOS_P(CellValues[0,:])
	states[1,3] = eos.IDEAL_EOS_P(CellValues[1,:])

	#Store rho*e
	states[0,4] = CellValues[0,2]
	states[1,4] = CellValues[1,2]

	return states


def Compute_AUSMPlus_Cell_States(CellValues):
	"""###################################################################################
	# The purpose of this function is to compute a set of state variables for a cell to be
	# used in the HLLC flux calculation.
	#
	# Author: Christopher Neal
	#
	# Date: 06-18-2015
	# Updated: 06-18-2015
	#
	######################################################################################
	#
	# Information: This function takes in a 2x3 array of cell left & right cell conservative
	#              variables and outputs a 2x5 array that contains the necessary cell variables
	#              for the HLLC flux scheme to estimate a flux. A combination of primitive and
	#              conserved variables.
	"""

	#DEBUG
        if(Debug_Flag==1):
                print 'Left Density = ' + str(CellValues[0,0])
                print 'Left Momentum = ' + str(CellValues[0,1])
                print 'Left Energy = ' + str(CellValues[0,2])

                print 'Right Density = ' + str(CellValues[1,0])
                print 'Right Momentum = ' + str(CellValues[1,1])
                print 'Right Energy = ' + str(CellValues[1,2])


	#Create an array to hold left and right state data
	states = np.zeros((2,5))

	#CellValues is a 2x3 array containing the conserved variable information for the
	#left and right states CellValues(0,:)=>Left State & CellValues(1,:)=>Right State

	#States is a 2x5 array that contains the relevant state information for the left
	#and right states that is computed from the conserved cell values i.e. primitive values.
	#the first row is states from the left cell & the second row is right cell states

	#Store the velocities of the states
	states[0,0] = CellValues[0,1]/CellValues[0,0] #Left velocity
	states[1,0] = CellValues[1,1]/CellValues[1,0] #Right velocity

	#Compute speed of sound(SS)
	states[0,1] = sqrt( gamma*eos.IDEAL_EOS_P(CellValues[0,:])/CellValues[0,0] ) #Left SS
	states[1,1] = sqrt( gamma*eos.IDEAL_EOS_P(CellValues[1,:])/CellValues[1,0] ) #Right SS
	
	#Store density
	states[0,2] = CellValues[0,0]
	states[1,2] = CellValues[1,0]

	#Store Pressure
	states[0,3] = eos.IDEAL_EOS_P(CellValues[0,:])
	states[1,3] = eos.IDEAL_EOS_P(CellValues[1,:])

	#Store rho*e
	states[0,4] = CellValues[0,2]
	states[1,4] = CellValues[1,2]


	return states



def Compute_ROE_Cell_States(CellValues):
	"""###################################################################################
	# The purpose of this function is to compute a set of state variables for a cell to be
	# used in the ROE flux calculation.
	#
	# Author: Christopher Neal
	#
	# Date: 10-22-2015
	# Updated: 10-22-2015
	#
	######################################################################################
	#
	# Information: This function takes in a 2x3 array of cell left & right cell conservative
	#              variables and outputs a 2x3 array that contains the necessary cell variables
	#              for the ROE flux scheme to estimate a flux. A combination of primitive and
	#              conserved variables.
	"""

	#DEBUG
        if(Debug_Flag==1):
                print 'Left Density = ' + str(CellValues[0,0])
                print 'Left Momentum = ' + str(CellValues[0,1])
                print 'Left Energy = ' + str(CellValues[0,2])

                print 'Right Density = ' + str(CellValues[1,0])
                print 'Right Momentum = ' + str(CellValues[1,1])
                print 'Right Energy = ' + str(CellValues[1,2])


	#Create an array to hold left and right state data
	states = np.zeros((2,3))

	#CellValues is a 2x3 array containing the conserved variable information for the
	#left and right states CellValues(0,:)=>Left State & CellValues(1,:)=>Right State

	#States is a 2x3 array that contains the relevant state information for the left
	#and right states that is computed from the conserved cell values i.e. primitive values.
	#the first row is states from the left cell & the second row is right cell states

	#Store the velocities of the states
	states[0,0] = CellValues[0,1]/CellValues[0,0] #Left velocity
	states[1,0] = CellValues[1,1]/CellValues[1,0] #Right velocity

	#Store density
	states[0,1] = CellValues[0,0]
	states[1,1] = CellValues[1,0]

	#Store h(enthalpy)
	states[0,2] = ( CellValues[0,2] + eos.IDEAL_EOS_P(CellValues[0,:])  )/CellValues[0,0]
	states[1,2] = ( CellValues[1,2] + eos.IDEAL_EOS_P(CellValues[1,:])  )/CellValues[1,0]


	return states

