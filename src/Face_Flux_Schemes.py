#######################################################################################
# This module contains the implementations of the different flux schemes that are used
# by the code.
#
# Author: Christopher Neal
#
# Date: 06-18-2015
# Updated: 06-18-2015
#
#######################################################################################

def hllc_flux(CellValues):
	"""################################################################################
	# The purpose of this function is to compute an estimate for the flux across a cell
	# face using the HLLC flux approximation method.
	#
	# Author: Christopher Neal
	#
	# Date: 06-18-2015
	# Updated: 06-18-2015
	#
	###################################################################################
	#
	# Information: This function takes in the states of the left and right cells and
	#              outputs the HLLC flux estimate for the face.
	"""
	from math import sqrt
	from Input_Params import *
	import numpy as np
	import Compute_States as CS
	import Compute_Cell_Flux as CF

	#The CellValues array is a 2x3 array containing the cell information for the
	#left and right adjacent cells. The cell information is the conserved variables
	
	#NOTE: Lowercase u variables correspond to x direction velocity. Uppercase U 
	#variables correspond to a 3 element vector containing conserved variables or
	#some quanity that is related to them.

	# Compute variables that are needed for the flux calculation and store in array
	states = CS.Compute_HLLC_Cell_States(CellValues)
	
	#Store states into named variables to make calculations easier to read
	uL = states[0,0]
	uR = states[1,0]

	aL = states[0,1]
	aR = states[1,1]

	rhoL = states[0,2]
	rhoR = states[1,2]

	pL = states[0,3]
	pR = states[1,3]

	EL = states[0,4] #Note E is simply rho*e. In the literature for HLLC they use E instead.
	ER = states[1,4]

	#DEBUG
	if(Debug_Flag == 1):
		print '\nLeft & Right States Cell Face'
		print 'uL = ' + str(uL)
		print 'uR = ' + str(uR)
		print 'aL = ' + str(aL)
		print 'aR = ' + str(aR)
		print 'rhoL = ' + str(rhoL)
		print 'rhoR = ' + str(rhoR)
		print 'pL = ' + str(pL)
		print 'pR = ' + str(pR)
	


	if(Wave_Estimate == 0): #Use TORO 1994 Estimate (Slide 31)

		#Compute average quantities between left and right states
		rho_bar = 0.5*(rhoL+rhoR)
		a_bar = 0.5*(aL+aR)

		Ppvrs = 0.5*(pL+pR) - 0.5*(uR-uL)*rho_bar*a_bar

		#Compute pressure estimate for middle pressure from equation
		p_star = max(0,Ppvrs)

		if(p_star <= pL):
			qL = 1.0	
		else:
			qL = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pL - 1))		

		SL = uL - aL*qL


		if(p_star <= pR):
                	qR = 1.0
        	else:
                	qR = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pR - 1))

        	SR = uR + aR*qR

	elif(Wave_Estimate == 1): #Use Min-Mod Selection

		SL = min(uL-aL,uR-aR)
		SR = max(uL+aL,uR+aR)

	#Compute middle wave speed
	S_star = (pR-pL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR))/(rhoL*(SL-uL) - rhoR*(SR-uR))

	#DEBUG
	if(Debug_Flag == 1):
		print 'SL = ' + str(SL)
		print 'SR = ' + str(SR)
		print 'S* = ' + str(S_star)	


	if(SL>=0):

		#The flux is just the value of the left side flux
		FluxEstimates = CF.Compute_Cell_Fluxes(CellValues[0,:])

		return FluxEstimates

	elif ( SL <=0 and S_star >=0):

		#Flux is given by the value of F_Star_L
		U_Star_L = np.zeros(3)

		U_Star_L[0] = 1
		U_Star_L[1] = S_star
		U_Star_L[2] = (EL/rhoL) + (S_star-uL)*(S_star + (pL)/(rhoL*(SL-uL)) )

		U_Star_L = ( rhoL*(SL-uL)/(SL-S_star) )*U_Star_L

		FL = CF.Compute_Cell_Fluxes(CellValues[0,:])

		F_Star_L = FL + SL*(U_Star_L - CellValues[0,:])

		FluxEstimates = F_Star_L

		return FluxEstimates		


	elif (S_star <= 0 and SR >=0 ):
	
		#Flux is given by the value of F_Star_R
                U_Star_R = np.zeros(3)

                U_Star_R[0] = 1
                U_Star_R[1] = S_star
                U_Star_R[2] = (ER/rhoR) + (S_star-uR)*(S_star + (pR)/(rhoR*(SR-uR)) )

                U_Star_R = ( rhoR*(SR-uR)/(SR-S_star) )*U_Star_R

                FR = CF.Compute_Cell_Fluxes(CellValues[1,:])

                F_Star_R = FR + SR*(U_Star_R - CellValues[1,:])

                FluxEstimates = F_Star_R
		
		return FluxEstimates

	elif (SR <= 0 ):

		#The flux is just the value of the right side flux
                FluxEstimates = CF.Compute_Cell_Fluxes(CellValues[1,:])
		
		return FluxEstimates
		
	else:
		print 'Error in Flux calculations in function hllc_flux'



def ausmPlus_flux(CellValues):
	"""################################################################################
	# The purpose of this function is to compute an estimate for the flux across a cell
	# face using the AUSM+ flux approximation method.
	#
	# Author: Christopher Neal
	#
	# Date: 09-30-2015
	# Updated: 09-30-2015
	#
	###################################################################################
	#
	# Information: This function takes in the states of the left and right cells and
	#              outputs the AUSM+ flux estimate for the face.
	"""
	from math import sqrt
	from Input_Params import *
	import numpy as np
	import Compute_States as CS
	import Compute_Cell_Flux as CF

	#The CellValues array is a 2x3 array containing the cell information for the
	#left and right adjacent cells. The cell information is the conserved variables
	
	#NOTE: Lowercase u variables correspond to x direction velocity. Uppercase U 
	#variables correspond to a 3 element vector containing conserved variables or
	#some quanity that is related to them.

	# Compute variables that are needed for the flux calculation and store in array
	states = CS.Compute_ausmPlus_Cell_States(CellValues)
	
	#Store states into named variables to make calculations easier to read
	uL = states[0,0]
	uR = states[1,0]

	aL = states[0,1]
	aR = states[1,1]

	rhoL = states[0,2]
	rhoR = states[1,2]

	pL = states[0,3]
	pR = states[1,3]

	EL = states[0,4] #Note E is simply rho*e. In the literature for HLLC they use E instead.
	ER = states[1,4]

	#DEBUG
	if(Debug_Flag == 1):
		print '\nLeft & Right States Cell Face'
		print 'uL = ' + str(uL)
		print 'uR = ' + str(uR)
		print 'aL = ' + str(aL)
		print 'aR = ' + str(aR)
		print 'rhoL = ' + str(rhoL)
		print 'rhoR = ' + str(rhoR)
		print 'pL = ' + str(pL)
		print 'pR = ' + str(pR)
	


	if(Wave_Estimate == 0): #Use TORO 1994 Estimate (Slide 31)

		#Compute average quantities between left and right states
		rho_bar = 0.5*(rhoL+rhoR)
		a_bar = 0.5*(aL+aR)

		Ppvrs = 0.5*(pL+pR) - 0.5*(uR-uL)*rho_bar*a_bar

		#Compute pressure estimate for middle pressure from equation
		p_star = max(0,Ppvrs)

		if(p_star <= pL):
			qL = 1.0	
		else:
			qL = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pL - 1))		

		SL = uL - aL*qL


		if(p_star <= pR):
                	qR = 1.0
        	else:
                	qR = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pR - 1))

        	SR = uR + aR*qR

	elif(Wave_Estimate == 1): #Use Min-Mod Selection

		SL = min(uL-aL,uR-aR)
		SR = max(uL+aL,uR+aR)

	#Compute middle wave speed
	S_star = (pR-pL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR))/(rhoL*(SL-uL) - rhoR*(SR-uR))

	#DEBUG
	if(Debug_Flag == 1):
		print 'SL = ' + str(SL)
		print 'SR = ' + str(SR)
		print 'S* = ' + str(S_star)	


	if(SL>=0):

		#The flux is just the value of the left side flux
		FluxEstimates = CF.Compute_Cell_Fluxes(CellValues[0,:])

		return FluxEstimates

	elif ( SL <=0 and S_star >=0):

		#Flux is given by the value of F_Star_L
		U_Star_L = np.zeros(3)

		U_Star_L[0] = 1
		U_Star_L[1] = S_star
		U_Star_L[2] = (EL/rhoL) + (S_star-uL)*(S_star + (pL)/(rhoL*(SL-uL)) )

		U_Star_L = ( rhoL*(SL-uL)/(SL-S_star) )*U_Star_L

		FL = CF.Compute_Cell_Fluxes(CellValues[0,:])

		F_Star_L = FL + SL*(U_Star_L - CellValues[0,:])

		FluxEstimates = F_Star_L

		return FluxEstimates		


	elif (S_star <= 0 and SR >=0 ):
	
		#Flux is given by the value of F_Star_R
                U_Star_R = np.zeros(3)

                U_Star_R[0] = 1
                U_Star_R[1] = S_star
                U_Star_R[2] = (ER/rhoR) + (S_star-uR)*(S_star + (pR)/(rhoR*(SR-uR)) )

                U_Star_R = ( rhoR*(SR-uR)/(SR-S_star) )*U_Star_R

                FR = CF.Compute_Cell_Fluxes(CellValues[1,:])

                F_Star_R = FR + SR*(U_Star_R - CellValues[1,:])

                FluxEstimates = F_Star_R
		
		return FluxEstimates

	elif (SR <= 0 ):

		#The flux is just the value of the right side flux
                FluxEstimates = CF.Compute_Cell_Fluxes(CellValues[1,:])
		
		return FluxEstimates
		
	else:
		print 'Error in Flux calculations in function hllc_flux'




def roe_flux(CellValues):
	"""################################################################################
	# The purpose of this function is to compute an estimate for the flux across a cell
	# face using the ROE flux approximation method.
	#
	# Author: Christopher Neal
	#
	# Date: 10-22-2015
	# Updated: 10-22-2015
	#
	###################################################################################
	#
	# Information: This function takes in the states of the left and right cells and
	#              outputs the ROE flux estimate for the face.
	"""
	from math import sqrt
	from Input_Params import *
	import numpy as np
	import Compute_States as CS
	import Compute_Cell_Flux as CF

	#The CellValues array is a 2x3 array containing the cell information for the
	#left and right adjacent cells. The cell information is the conserved variables
	
	#NOTE: Lowercase u variables correspond to x direction velocity. Uppercase U 
	#variables correspond to a 3 element vector containing conserved variables or
	#some quanity that is related to them.

	# Compute variables that are needed for the flux calculation and store in array
	states = CS.Compute_ROE_Cell_States(CellValues)
	
	#Store states into named variables to make calculations easier to read
	uL = states[0,0]
	uR = states[1,0]

	rhoL = states[0,1]
	rhoR = states[1,1]

	hL = states[0,2] #Enthalpy.
	hR = states[1,2]

	#DEBUG
	if(Debug_Flag == 1):
		print '\nLeft & Right States Cell Face'
		print 'uL = ' + str(uL)
		print 'uR = ' + str(uR)
		print 'rhoL = ' + str(rhoL)
		print 'rhoR = ' + str(rhoR)
		print 'hL = ' + str(hL)
		print 'hR = ' + str(hR)
	

	#Compute density averaged quantities
	hBar = (sqrt(rhoL)*hL + sqrt(rhoR)*hR)/(sqrt(rhoL) + sqrt(rhoR))
	uBar = (sqrt(rhoL)*uL + sqrt(rhoR)*uR)/(sqrt(rhoL) + sqrt(rhoR))
	cBar = sqrt((gamma-1)*(hBar - 0.5*uBar*uBar))

	#Define characteristic speeds
	lambda1 = uBar - cBar
	lambda2 = uBar
	lambda3 = uBar + cBar	

	#Create r_p vectors
	r_1 = np.zeros(3)
	r_1[0] = 1
	r_1[1] = uBar - cBar
	r_1[2] = hBar - uBar*cBar

	r_2 = np.zeros(3)
	r_2[0] = 1
	r_2[1] = uBar
	r_2[2] = 0.5*uBar*uBar

	r_3 = np.zeros(3)
        r_3[0] = 1
        r_3[1] = uBar + cBar
        r_3[2] = hBar + uBar*cBar

	#Compute the r^p vectors for computing alpha_p
	r1 = np.zeros(3)
	r1[0] = (uBar/(4*cBar))*(2 + (gamma-1)*(uBar/cBar))
	r1[1] = (-1.0/(2.0*cBar))*(1 + (gamma-1)*(uBar/cBar))
	r1[2] = 0.5*(gamma-1)*(1/(cBar*cBar))	
	
	r2 = np.zeros(3)
        r2[0] = 1 - 0.5*(gamma-1)*((uBar*uBar)/(cBar*cBar))
        r2[1] = (gamma-1)*(uBar/(cBar*cBar))
        r2[2] = -(gamma-1)*(1.0/(cBar*cBar))

	r3 = np.zeros(3)
        r3[0] = -(uBar/(4*cBar))*(2 - (gamma-1)*(uBar/cBar))
        r3[1] = (1.0/(2.0*cBar))*(1 - (gamma-1)*(uBar/cBar))
        r3[2] = 0.5*(gamma-1)*(1/(cBar*cBar))

	alpha1 = np.inner(r1 , CellValues[1,:] - CellValues[0,:])
	alpha2 = np.inner(r2 , CellValues[1,:] - CellValues[0,:])
	alpha3 = np.inner(r3 , CellValues[1,:] - CellValues[0,:])


	Term1 = abs(lambda1)*alpha1*r_1
	Term2 = abs(lambda2)*alpha2*r_2
	Term3 = abs(lambda3)*alpha3*r_3
	FluxAverage = 0.5*(CF.Compute_Cell_Fluxes(CellValues[1,:]) + CF.Compute_Cell_Fluxes(CellValues[0,:]))	

	FluxEstimates = FluxAverage -0.5*(Term1 + Term2 + Term3)
	return FluxEstimates		


