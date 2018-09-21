"""
This module contains functions for computing important cell state values for different
flux schems.
"""
from math import sqrt
import numpy as np

class CellStates(object):
	"""
	The purpose of this class is to compute a set of state variables for a cell to be
	used in flux calculations.
	"""
    def __init__(self):
        pass

    def compute_states(self, adjacent_cell_values):
            raise NotImplementedError

    def print_cell_values(self, adjacent_cell_values):
        print 'Left Density = ' + str(adjacent_cell_values['left_cell']['rho'])
        print 'Left Momentum = ' + str(adjacent_cell_values['left_cell']['rho*u'])
        print 'Left Energy = ' + str(adjacent_cell_values['left_cell']['rho*e'])

        print 'Right Density = ' + str(adjacent_cell_values['right_cell']['rho'])
        print 'Right Momentum = ' + str(adjacent_cell_values['right_cell']['rho*u'])
        print 'Right Energy = ' + str(adjacent_cell_values['right_cell']['rho*e'])


class HLLCStates(CellStates):
    def __init__(self):
        pass

    def compute_states(self, adjacent_cell_values):
	#States is a dictionary that contains the relevant state information for the left
	#and right states that is computed from the conserved cell values i.e. primitive values.
    states = {}

	#Store the velocities of the states
	states['u_left'] = adjacent_cell_values['left_cell']['rho*u'] / adjacent_cell_values['left_cell']['rho'] #Left velocity
	states['u_right'] =  adjacent_cell_values['right_cell']['rho*u'] / adjacent_cell_values['right_cell']['rho'] #Right velocity

	#Compute speed of sound(SS)
	states['soundspeed_left'] = self.eos.soundspeed(adjacent_cell_values['left_cell']) #Left SS
	states['soundspeed_right'] = self.eos.soundspeed(adjacent_cell_values['right_cell']) #Right SS
	
	#Store density
	states['rho_left'] = adjacent_cell_values['left_cell']['rho']
	states['rho_right'] = adjacent_cell_values['right_cell']['rho']

	#Store Pressure
	states['pressure_left'] = self.eos.pressure(adjacent_cell_values['left_cell'])
	states['pressure_right'] = self.eos.pressure(adjacent_cell_values['right_cell'])

	#Store rho*e
	states['rho*e_left']= adjacent_cell_values['left_cell']['rho*e']
	states['rho*e_right']= adjacent_cell_values['right_cell']['rho*e']

	return states


class ROEStates(CellStates):
    def __init__(self):
        pass

    def compute_cell_states(self, adjacent_cell_values):
        """
        The purpose of this function is to compute a set of state variables for a cell to be
        used in the ROE flux calculation.
        """
        #States is a dictionary that contains the relevant state information for the left
        #and right states that is computed from the conserved cell values i.e. primitive values.
        states = {}

        #Store the velocities of the states
        states['u_left'] = adjacent_cell_values['left_cell']['rho*u'] / adjacent_cell_values['left_cell']['rho'] #Left velocity
        states['u_right'] = adjacent_cell_values['right_cell']['rho*u'] / adjacent_cell_values['right_cell']['rho'] #Right velocity

        #Store density
        states['rho_left'] = adjacent_cell_values['left_cell']['rho']
        states['rho_right'] = adjacent_cell_values['right_cell']['rho']

        #Store h(enthalpy)
        states['enthalpy_left'] = ( adjacent_cell_values['left_cell']['rho*e'] + self.eos.pressure(adjacent_cell_values['left_cell'])  ) / adjacent_cell_values['left_cell']['rho']
        states['enthalpy_right'] = ( adjacent_cell_values['right_cell']['rho*e'] + self.eos.pressure(adjacent_cell_values['right_cell'])  ) / adjacent_cell_values['right_cell']['rho']


        return states

