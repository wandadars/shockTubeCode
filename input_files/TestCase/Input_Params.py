#######################################################################################
# The purpose of this Python module is to store all of the global parameters that will
# be frequenty used by the program.
#
#
# Author: Christopher Neal
#
# Date:         06-18-2015
# Updated:      06-23-2015
#
#######################################################################################

#Geometric Parameters for Problem
x0 = 0.4        	# Initial x coordinate of domain
xL = 0.6        	# Max x coordinate of domain
Nx = 200		# Number of cells to create in the domain
dt = 2.5e-4		# Timestep
Nt = 1000       	# Number of timesteps to take

#Solution Printing
Print_Time = 2.5e-2

#Solution Plotting real-time
Real_Time_Plot = True


#Gas Properties
gamma = 1.4		#Specific Heat Ratio
R_gas = 287.0 		# J/kgK

#File Output Format
Output_Format = 1			# 0 - regular data format, 1 - VTK Legacy
VTK_Grid_Vert_Cells = 200		# Number of duplicate cells to use for VTK output


#Numerics
Wave_Estimate = 1	#0-Method 1, 1-TORO's Method


#Initialization Parameters
Init_Type = 0
X_Loc_1 = 0.5

#Left Side
u1 = 0.0
P1 = 1.0
rho1 = 1.0

#Right Side
u2 = 0.0
P2 = 0.1
rho2 = 0.125


#Boundary Conditions
uLeft = 0.0       	# meters/second
PLeft = 1.0	 	# Pascals
TLeft = 0.0034843     	# Kelvin

uRight = 0.0       	# meters/second
PRight = 0.1	  	# Pascals
TRight = 0.027875     	# Kelvin


#Debug Section( 0=> no debugging, 1=> debugging statements on )
Debug_Flag = 0


#Compute dx from inputs
dx = (xL-x0)/float(Nx)


