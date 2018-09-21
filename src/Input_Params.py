#geometric parameters for problem
x0  0.0     # Initial x coordinate of domain
xL  1.0     # Max x coordinate of domain
Nx  800	    # Number of cells to create in the domain
dt  2.5e-4  # Timestep
Nt  2000    # Number of timesteps to take

#solution printing frequency
Print_Time  2.5e-3

#solution plotting real-time
Real_Time_Plot False

#Gas Properties
eos IDEAL
gamma  1.4		#Specific Heat Ratio
R_gas  287.0 		# J/kgK

#File Output Format
output_format  1			# 0 - regular data format, 1 - VTK Legacy
VTK_Grid_Vert_Cells  100   # Number of vertical cells to duplicate solution onto 

#Numerics
wave_estimate  1	#0-Method 1, 1-TORO's Method
flux_scheme HLLC	#0-HLLC, 1-ROE, 2-AUSM+, 3-WENO5

#Initialization Parameters
Init_Type  0
X_Loc_1  0.5

#Left Side
u1  0.0
P1  1.0
rho1  1.0

#Right Side
u2  0.0
P2  0.1
rho2  0.125

#Boundary Conditions
uLeft  0.0       	# meters/second
PLeft  1.0	 	# Pascals
TLeft  0.0034843     	# Kelvin

uRight  0.0       	# meters/second
PRight  0.1	  	# Pascals
TRight  0.0027875     	# Kelvin


#Debug Section( 0=> no debugging, 1=> debugging statements on )
debug_flag  0



