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
import numpy as np
import matplotlib.pyplot as plt
import time

from pylab import *

from Input_Params import *
from Initialize import *
from RK_Method import *
from WriteSolution import *
from Compute_Primitives import *



class Main(object):
    def __init__(self, input_file_name):
        self.initial_program_state = {}
        input_parser = ip.InputFileParser(input_file_name)
        self.initial_program_state['input_parser'] = input_parser

        if input_parser.user_input_data['flux_type'].lower() == 'roe':
            flux_scheme = fs.RoeFluxScheme()



#Create 1D coordinate vector and initialize for plotting solution
x = np.zeros(Nx)
x[0] = dx/2.0
for i in range(1, Nx, 1):
	x[i] = x[i-1] + dx

# Create array for storing the conserved solution variables in the domain. 
# We are solving a system of 3 equations, so 3 columns are allocated.
G = np.zeros((Nx,3)) # Array for holding the conserved cell variables(rho, rho*u, rho*E)


# Initialize the domain
G = Initialize_Domain()

#Write Initial Domain Data
write_solution(x, G, 0, 0*dt)


if(Real_Time_Plot == True):
	#Plotting Setup
	plt.ion()
	fig = plt.figure( 1,figsize = (10,10) )

#Initialize printing variable
Elapsed_Time = 0.0

# Time Stepping Loop
for i in range(0,Nt,1):
	G_New = RK_3_Stepping(G)	

	if(Real_Time_Plot == True):
		Prim = np.zeros((Nx,4))
		Prim = compute_primitives(G_New)

		#Plot Solution
		if(i == 0 ):
			#Set up the plots
			ax = fig.add_subplot(221)

			#Plot Density
			line1, = ax.plot(x,Prim[:,0],'r-')
			plt.xlabel('Position, x (meters)')
			plt.ylabel('Density, rho (kg/m^3)')
			plt.axis((0.0,1.0,0.1,1.05))

			#Plot Velocity
			ax2 = fig.add_subplot(222)
			line2, = ax2.plot(x,Prim[:,1],'r-')
                        plt.xlabel('Position, x (meters)')
                        plt.ylabel('Velocity, u (m/s)')
                        plt.axis((0.0,1.0,0.0,1.0))

			#Plot Pressure
			ax3 = fig.add_subplot(223)
                        line3, = ax3.plot(x,Prim[:,2],'r-')
                        plt.xlabel('Position, x (meters)')
                        plt.ylabel('Pressure, P (Pascal )')
                        plt.axis((0.0,1.0,0.05,1.05))


			#Plot Temperature
                        ax4 = fig.add_subplot(224)
                        line4, = ax4.plot(x,Prim[:,3],'r-')
                        plt.xlabel('Position, x (meters)')
                        plt.ylabel('Temperature, T (Kelvin )')
                        plt.axis((0.0,1.0,0.0,0.005))


		else:

			line1.set_ydata(Prim[:,0])
			line2.set_ydata(Prim[:,1])
			line3.set_ydata(Prim[:,2])
			line4.set_ydata(Prim[:,3])

			fig.canvas.draw()
			time.sleep(0.001)
	


	#Update conserved variables i.e. make a step in time
	G = G_New

	#Increment Elapsed Time
	Elapsed_Time = Elapsed_Time + dt

	if( Elapsed_Time >= Print_Time ):
		write_solution(x,G,i+1,(i+1)*dt)
		Elapsed_Time = 0.0 #reset elapsed time
		
	print("Completed Timestep:\t%d\n"%(i+1))


