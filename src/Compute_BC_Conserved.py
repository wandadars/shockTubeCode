#######################################################################################
# The purpose of this function is to compute an estimate for the value of the 
# conserved variable on the boundary.
#
# Author: Christopher Neal
#
# Date: 06-18-2015
# Updated: 06-18-2015
#
#######################################################################################


def compute_BC_conserved(BC_Num):

	from Input_Params import *
	import numpy as np

	BC_Flux = np.zeros(3)	

	if(BC_Num == 0): #Left boundary conserved variables

		BC_Flux[0] = (PLeft/(R_gas*TLeft))	#rho
		BC_Flux[1] = BC_Flux[0]*uLeft		#rho*u
		#rho*e
		BC_Flux[2] = BC_Flux[0]*((R_gas/(gamma-1))*TLeft + 0.5*uLeft*uLeft)
	
	elif(BC_Num == 1): #Right boundary conserved variables
	
		BC_Flux[0] = (PRight/(R_gas*TRight))     #rho
                BC_Flux[1] = BC_Flux[0]*uRight            #rho*u

		#rho*e
                BC_Flux[2] = BC_Flux[0]*((R_gas/(gamma-1))*TRight + 0.5*uRight*uRight)


	else:
		print 'Error: Improper BC Number specified for compute_BC_conserved() function'

	
	#DEBUG
	if(Debug_Flag == 1):
		
		if(BC_Num == 0):
			print '\nBC Conserved for Left Boundary: \n'
			print 'rho = ' + str(BC_Flux[0])
			print 'rho*u = ' + str(BC_Flux[1])
			print 'rho*e = ' + str(BC_Flux[2])

		elif(BC_Num == 1):
			print '\nBC Conserved for Right Boundary: \n'
                        print 'rho = ' + str(BC_Flux[0])
                        print 'rho*u = ' + str(BC_Flux[1])
                        print 'rho*e = ' + str(BC_Flux[2])


	return BC_Flux

