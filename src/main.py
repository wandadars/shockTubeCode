"""
The purpose of this Python script is to generate the solution to the standard 1D SOD
shock tube problem. It solves the unsteady compressible Euler equations.
 Information:
 This particular solver uses the following schemes:
		1.) 3 Stage Runge-Kutta time discretization scheme.
		2.) HLLC flux approximation for the cell face fluxes.
		3.) Finite volume treatment of the governing equations.
"""

import sys

#Modules
import input_parser as ip
import flux_schemes as fs
import solver as slv
import equation_of_state as eos
import rk_method

class Main(object):
    def __init__(self, input_file_name):
        self.initial_program_state = {}
        input_parser = ip.InputFileParser(input_file_name)
        self.initial_program_state['input_parser'] = input_parser

        if input_parser.user_input_data['flux_type'].lower() == 'roe':
            flux_calculator = fs.ROEFluxScheme()
        if input_parser.user_input_data['flux_type'].lower() == 'hllc':
            flux_calculator = fs.HLLCFluxScheme()

        if input_parser.user_input_data['eos'].lower() == 'ideal':
            r_gas = float(input_parser.user_input_data['gamma'])
            gamma = float(input_parser.user_input_data['R_gas'])
            eos_object = eos.IdealEOS(gamma, r_gas)

        if input_parser.user_input_data['init_type'].lower() == '0':
            flow_initializer = slv.ShockTubeInitializer(input_parser, eos_object)

        if input_parser.user_input_data['time_integrator'].lower() == 'rk3':
            time_integrator = rk_method.RK3Integrator(flux_scheme)

        solution = slv.Solution(input_parser, flow_initializer, time_integrator, solution_writer_object)
        
        self.solver = slv.Solver(input_parser, solution)

    def run(self):
        self.solver.begin_timestepping()

file_input = sys.argv[1]
main = Main(file_input)
Main.run()

