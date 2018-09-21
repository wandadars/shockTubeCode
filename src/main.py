"""
The purpose of this Python script is to generate the solution to the standard 1D SOD
shock tube problem. It solves the unsteady compressible Euler equations.
 Information:
 This particular solver uses the following schemes:
		1.) 3 Stage Runge-Kutta time discretization scheme.
		2.) HLLC flux approximation for the cell face fluxes.
		3.) Finite volume treatment of the governing equations.
"""

#Modules
import input_parser as ip
import flux_schemes as fs


class Main(object):
    def __init__(self, input_file_name):
        self.initial_program_state = {}
        input_parser = ip.InputFileParser(input_file_name)
        self.initial_program_state['input_parser'] = input_parser

        if input_parser.user_input_data['flux_type'].lower() == 'roe':
            flux_scheme = fs.ROEFluxScheme()
        if input_parser.user_input_data['flux_type'].lower() == 'hllc':
            flux_scheme = fs.HLLCFluxScheme()





main = Main()

Main.run()

