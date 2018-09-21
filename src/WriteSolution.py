from Input_Params import *
from Compute_Primitives import *
import numpy as np
import os
from cmath import sqrt

def write_solution(x,CellValues,iteration,currentTime):
    """
    The purpose of this function is to output solution data to a file for post processing
    at a later time.
        Information: This function takes in the conserved variables array, G
	"""
	# CellValues is a Nx X 3 element vector that contains the conserved variables.

	# Before writing solution check to see if the output directory exists
	path = os.getcwd()
	if(os.path.isdir(path+'/output') == False): #Directory does not exist. Create it
		os.makedirs('output')
	
	OutPath = path+'/output'
	
	#Enter into output directory and print solution
	os.chdir(OutPath)

	if(Output_Format == 0):
	  #Create filename
	  filename = 'solution_'+ str(iteration) +'_'+ str(currentTime)

	  #Open output file
	  target = open(filename, 'w')

	  #Compute primitives for output
	  Prim = np.zeros((Nx,4)) #rho, u,P,T
	  Prim = compute_primitives(CellValues)

	  #Output to file
	  for i in range(0,Nx,1):
	    outputString = str(x[i])+ '\t' +str(Prim[i,0])+'\t'+str(Prim[i,1])+'\t'+str(Prim[i,2])+'\t'+str(Prim[i,3])+'\n'
	    target.write(outputString)

	  target.close()


	elif(Output_Format == 1):	#VTK Legacy format
	  #Determine appropriate number of padding zeros for output
	  pad = '00000'
	  if(iteration < 10):
	  	pad = '0000'
	  elif(iteration < 100):
	  	pad = '000'
	  elif(iteration < 1000):
	  	pad = '00'
	  elif(iteration < 10000):
	  	pad = '0'
	  else:
	  	pad = ''

	  #Create filename
      filename = 'solution_'+ pad + str(iteration)+'.vtk'

      #Open output file
      target = open(filename, 'w')

	  #Write the first line of the VTK file
	  target.write('# vtk DataFile Version 3.1 \n')
	  target.write('Output for 1D Compressible Euler Equation Solver - Chris Neal \n')
	
	  #Specify output format
	  target.write('ASCII \n')
	
	  #Specify the type of data set
	  target.write('DATASET UNSTRUCTURED_GRID \n\n')
	  
	  #Write the number of points
	  tempString = 'POINTS ' + str(Nx*VTK_Grid_Vert_Cells) + ' FLOAT\n'
	  target.write(tempString)

	  #Create the x and y arrays that hold coordinates
	  x=np.zeros(Nx)
	  y=np.zeros(VTK_Grid_Vert_Cells)
	  y[0] = dx/2.0
	  x[0] = dx/2.0
	  for i in range(1, Nx,1):
          	x[i] = x[i-1] + dx

	  for i in range(1,VTK_Grid_Vert_Cells,1):
	  	y[i] = y[i-1] + dx

	  #Write the point coordinates
	  for j in range(0,VTK_Grid_Vert_Cells,1):
	  	for i in range(0,Nx,1):
	  		tempString = str(x[i])+ '\t' + str(y[j])+ '\t' + '0\n'
	  		target.write(tempString)

	  target.write('\n')

	  #Specify Cells Header
	  NumCells = (Nx-1)*(VTK_Grid_Vert_Cells-1)
	  tempString = 'CELLS ' + str( NumCells) + ' ' + str(NumCells*5)+'\n'
	  target.write(tempString)
	
	  #Write Cell Connectivity
	  for j in range(1,VTK_Grid_Vert_Cells,1):
	  	for i in range(1,Nx,1):

			point1 = str( i + (j-1)*Nx - 1 )
			point2 = str( i + (j-1)*Nx )
			point3 = str( i + j*Nx )
			point4 = str( i + j*Nx - 1 )

	  		tempString = '4\t'+point1+'\t'+point2+'\t'+point3+'\t'+point4+'\n' 
	  		target.write(tempString)

	  target.write('\n')

	  #Write Cell Type Data
	  tempString = 'CELL_TYPES ' + str(NumCells) + '\n'
	  target.write(tempString)
	  for i in range(0,NumCells,1):
	  	target.write('9\n')

	  target.write('\n')


	  #Compute Data to Ouput
	  #Compute primitives for output
          Prim = np.zeros((Nx,4)) #rho, u, P, T
          Prim = compute_primitives(CellValues)

	  #Write Point Data
	  tempString = 'POINT_DATA ' + str(Nx*VTK_Grid_Vert_Cells)+'\n'
	  target.write(tempString)

	  #WRITE DENSITY
	  target.write('SCALARS Density FLOAT 1 \n')
	  target.write('LOOKUP_TABLE default \n')
	 
	  for j in range(0,VTK_Grid_Vert_Cells,1):
	  	for i in range(0,Nx,1):
			tempString = str(Prim[i,0]) + '\n'
			target.write(tempString)

	  #WRITE VELOCITY
	  target.write('\nSCALARS Velocity FLOAT 1 \n')
          target.write('LOOKUP_TABLE default \n')

          for j in range(0,VTK_Grid_Vert_Cells,1):
                for i in range(0,Nx,1):
                        tempString = str(Prim[i,1]) + '\n'
                        target.write(tempString)

	  #WRITE PRESSURE
	  target.write('\nSCALARS Pressure FLOAT 1 \n')
          target.write('LOOKUP_TABLE default \n')

          for j in range(0,VTK_Grid_Vert_Cells,1):
                for i in range(0,Nx,1):
                        tempString = str(Prim[i,2]) + '\n'
                        target.write(tempString)

	  #WRITE TEMPERATURE
	  target.write('\nSCALARS Temperature FLOAT 1 \n')
          target.write('LOOKUP_TABLE default \n')

          for j in range(0,VTK_Grid_Vert_Cells,1):
                for i in range(0,Nx,1):
                        tempString = str(Prim[i,3]) + '\n'
                        target.write(tempString)

	  #WRITE SOUNDSPEED
          target.write('\nSCALARS SOUNDSPEED FLOAT 1 \n')
          target.write('LOOKUP_TABLE default \n')

          for j in range(0,VTK_Grid_Vert_Cells,1):
                for i in range(0,Nx,1):
			soundspeed = sqrt(gamma*R_gas*Prim[i,3])
                        tempString = str(soundspeed) + '\n'
                        target.write(tempString)

	  #WRITE SOUNDSPEED
          target.write('\nSCALARS MACH FLOAT 1 \n')
          target.write('LOOKUP_TABLE default \n')

          for j in range(0,VTK_Grid_Vert_Cells,1):
                for i in range(0,Nx,1):
                        mach = Prim[i,1]/sqrt(gamma*R_gas*Prim[i,3])
                        tempString = str(mach) + '\n'
                        target.write(tempString)



	  target.close()

	#Return back up into the working directory
	os.chdir(path)




