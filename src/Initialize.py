#######################################################################################
# The purpose of this Python script is to loop over the domain at time=0 and initialize
# all of the flow variables.             
#                                                                                      
#                                                                                      
# Author: Christopher Neal                                                             
#                                                                                      
# Date:         06-18-2015                                                             
# Updated:      06-18-2015                                                             
#                                                                                      
#######################################################################################
#                                                                                      

def Initialize_Domain():

	from Input_Params import *	
	import numpy as np
	
	G = np.zeros((Nx,3))
	
	#Create X Vector and initialize for plotting solution

	#DEBUG
	if(Debug_Flag==1):
		print '\n dx = '+str(dx)+ '\n'
		print 'Left Boundary Location: x0 = '+str(x0)+'\n'
		print 'Right Boundary Location: xL = '+str(xL)+'\n'

        x = np.zeros(Nx)
        x[0] = x0+dx/2.0
        for i in range(1, Nx,1):
		x[i] = x[i-1] + dx


	if(Init_Type == 0):

		#Neal Shock Initialization

		#Debug
		if(Debug_Flag == 1):
			print 'Left and Right Initializations are:\n'
			
			print 'Left Properties:\n'
			print 'rho1 = '+str(rho1)
			print 'rho1*u1 = '+str(rho1*u1)
			print 'rho1*e1 = '+str(rho1*( P1/(rho1*(gamma-1))) + 0.50*(u1*u1))		
		
			print 'Right Properties:\n'
			print 'rho2 = '+str(rho2)
                        print 'rho2*u2 = '+str(rho2*u2)
                        print 'rho2*e2 = '+str(rho2*( P2/(rho2*(gamma-1))) + 0.50*(u2*u2))

			print 'Shock is Initialized at:\n'
			print 'X_Shock = '+str(X_Loc_1)



		#Loop over all cells and assign initial conditions
		for i in range(0,Nx,1):

			if( x[i] <= X_Loc_1 ):

				G[i,0] = rho1
				G[i,1] = rho1*u1
				G[i,2] = rho1*( P1/(rho1*(gamma-1))) + 0.50*(u1*u1)
		
			else:
				G[i,0] = rho2
                        	G[i,1] = rho2*u2
                        	G[i,2] = rho2*( P2/(rho2*(gamma-1))) + 0.50*(u2*u2)

	#DEBUG
	if(Debug_Flag == 1): #Print Solution Vector
		print '\n Displaying Initialized Flowfield Vector\n'
		print '\t\t'+'X\t'+'rho'+'\t'+'rho*u'+'\t'+'rho*e\n'
		for i in range(0,Nx,1):
			print 'Cell '+str(i+1)+':\t\t'+str(x[i])+'\t'+ str(G[i,0]) + '\t' + str(G[i,1]) + '\t'+str(G[i,2]) + '\n'




	return G
	



