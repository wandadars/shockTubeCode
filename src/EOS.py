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

    def energy(cell_state):
        return cell_state['rho*e'] / cell_state['rho']

    def energy_PTU(pressure, temperature, velocity):
        return (self.r_gas * temperature / (self.gamma-1)) + 0.5*velocity**2
    
    def energy_PRU(pressure, density, velocity):
        return pressure / (density*(self.gamma-1)) + 0.5*velocity**2
    
    def pressure(cell_state):
	    # P = (gamma-1)*(rho*e - (1/2)*rho*u*u)
        return (self.gamma - 1.0) * (cell_state['rho*e'] - 0.5*cell_state['rho*u']**2 / cell_state['rho'])

    def density(cell_state)
        return cell_state['rho']

    def density_PT(pressure, temperature):
        return pressure / (self.r_gas * temperature)

    def soundspeed(cell_state):
        # a = sqrt(gamma*(gamma-1)*(e-0.5*u*u))
        return np.sqrt(self.gamma * (self.gamma-1) * (cell_state['rho*e']/cell_state['rho'] -0.5*cell_state['rho*u']**2 / cell_state['rho'])
