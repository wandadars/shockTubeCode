
class FaceFluxScheme(object):
    def __init__(self, eos):
        self.eos = eos

    def compute_flux(cell_values):
        raise NotImplementedError

    def compute_cell_fluxes(cell_value):
        """
        The purpose of this function is to compute the value of the flux variables at the cell
        centers using the conserved variables for the cell.

        Information: This function takes in a cell data 1x3 array that holds the 3 conserved
        variables.
            Output: A 1x3 array of the fluxes for mass, momentum,and energy
        """
        #Allocate the array to hold the flux data
        cell_flux = []

        #Store flux for conservation of mass equation i.e. rho*u
        cell_flux.append(cell_value['rho*u'])

        #Store flux for conservation of momentum equation i.e. rho*u^2 + p
        density = cell_value['rho']
        energy = cell_value['rho*e'] / cell_value['rho']
        velocity = cell_value['rho*u'] / cell_value['rho']
        cell_flux.append((density * velocity**2) + eos.pressure(density, energy, velocity))

        #Store flux for conservation of energy equation i.e. u*(rho*e + p)
        cell_flux.append(velocity * (density * energy + eos.pressure(density, energy, velocity)))

        return cell_flux

    
def compute_cell_face_fluxes(G):
    """
    The purpose of this function is to compute the values of the cell face fluxes the
    hllc method.

    Information: This function takes in a cell data Nx X 3 array that holds the 3 conserved
    variables for all cells in the domain.
        Output: A Nx X 3 X 2 array of the fluxes for mass, momentum,and energy at the west
        and east faces for each cell in the domain
    """

    #Allocate array to hold the conserved variables for the left and right cells that surround
    #the face.
    temp = np.zeros((2,3))

    #Compute fluxes at faces by looping over all cells and using a flux function

    #Create array to hold face flux values
    # the 2 is for the east and west fluxes for each cell.
    FaceFluxes = np.zeros((Nx,3,2)) # the 3 is because there are 3 conserved variables

    #The First Cell has a prescribed flux on its left face from the boundary
    #DEBUG
    if(Debug_Flag == 1):
        print 'Debug Info for Cell: '+ str(1) +'\n' 

    #Compute and store left face flux for the first cell
    temp[0,:] = BCONSERV.compute_BC_conserved(0)
    temp[1,:] = G[0,:]

    if(FluxSchemeInput == 0):       #HLLC Flux
        FaceFluxes[0,:,0] = FFS.hllc_flux(temp)

    elif(FluxSchemeInput == 1):     #ROE Flux
        FaceFluxes[0,:,0] = FFS.roe_flux(temp)


    #Compute and store right face flux for the first cell
    temp[0,:] = G[0,:]
    temp[1,:] = G[1,:]

    if(FluxSchemeInput == 0):   #HLLC Flux
        FaceFluxes[0,:,1] = FFS.hllc_flux(temp)

    elif(FluxSchemeInput == 1): #ROE Flux
        FaceFluxes[0,:,1] = FFS.roe_flux(temp)
    
    
    #Inner cells have 2 face fluxes that must be determined
    for j in range(1,Nx-1,1):
        #DEBUG
        if(Debug_Flag == 1):
            print 'Debug Info for Cell: '+ str(j+1) +'\n'

        #DEBUG
        if(Debug_Flag == 1):
            print 'Left Cell Face Flux:\n'

        #Compute and store left face flux
        if(FluxSchemeInput == 0):   #HLLC Flux
            #For this type of flux scheme the result of this will be the same as for the previous
            #cell's right face flux
            FaceFluxes[j,:,0] = FaceFluxes[j-1,:,1]
        elif(FluxSchemeInput == 1): #Roe Flux
            #For this type of flux scheme the result of this will be the same as for the previous
            #cell's right face flux
            FaceFluxes[j,:,0] = FaceFluxes[j-1,:,1]

        #DEBUG
        if(Debug_Flag == 1):
            print 'Right Cell Face Flux:\n'

        #Compute and store right face flux
        temp[0,:] = G[j,:]
        temp[1,:] = G[j+1,:]

        if(FluxSchemeInput == 0):   #HLLC Flux
            FaceFluxes[j,:,1] = FFS.hllc_flux(temp)
        elif(FluxSchemeInput == 1): #Roe Flux
            FaceFluxes[j,:,1] = FFS.roe_flux(temp)


    #The Last Cell has a prescribed flux on its right face from the boundary
    if(Debug_Flag == 1):
        print 'Debug Info for Cell: '+ str(Nx) +'\n'
        print 'Flux across Left Cell Face'
        

    #Compute and store left face flux
    temp[0,:] = G[Nx-2,:]   #Remember arrays go from 0 to Nx-1 for Nx elements
    temp[1,:] = G[Nx-1,:]

    if(FluxSchemeInput == 0):   #HLLC Flux
        FaceFluxes[Nx-1,:,0] = FFS.hllc_flux(temp)

    elif(FluxSchemeInput == 1): #Roe Flux
        FaceFluxes[Nx-1,:,0] = FFS.roe_flux(temp)

    #Compute and store right face flux
    temp[0,:] = G[Nx-1,:]
    temp[1,:] = BCONSERV.compute_BC_conserved(1)

    if(FluxSchemeInput == 0):       #HLLC Flux
        FaceFluxes[Nx-1,:,1] = FFS.hllc_flux(temp)

    elif(FluxSchemeInput == 1):     #Roe Flux
        FaceFluxes[Nx-1,:,1] = FFS.roe_flux(temp)


    return FaceFluxes



class HLLCFluxScheme(FaceFluxScheme):
    def __init__(self, eos


    def compute_flux(cell_values):
        """
        # The purpose of this function is to compute an estimate for the flux across a cell
        # face using the HLLC flux approximation method.
        #
        # Information: This function takes in the states of the left and right cells and
        #              outputs the HLLC flux estimate for the face.
        """
        #The CellValues array is a 2x3 array containing the cell information for the
        #left and right adjacent cells. The cell information is the conserved variables
        
        #NOTE: Lowercase u variables correspond to x direction velocity. Uppercase U 
        #variables correspond to a 3 element vector containing conserved variables or
        #some quanity that is related to them.

        # Compute variables that are needed for the flux calculation and store in array
        states = self.compute_cell_states(cell_states)
        
        #Store states into named variables to make calculations easier to read
        uL = states['u_left']
        uR = states['u_right']

        aL = states['soundspeed_left']
        aR = states['soundspeed_right']

        rhoL = states['rho_left']
        rhoR = states['rho_right']

        pL = states['pressure_left']
        pR = states['pressure_right']

        EL = states['rho*e_left'] #Note E is simply rho*e. In the literature for HLLC they use E instead.
        ER = states['rho*e_right']

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
        


        if Wave_Estimate == 0: #Use TORO 1994 Estimate (Slide 31)
            #Compute average quantities between left and right states
            rho_bar = 0.5 * (rhoL + rhoR)
            a_bar = 0.5*(aL + aR)
            Ppvrs = 0.5*(pL + pR) - 0.5 * (uR - uL) * rho_bar * a_bar

            #Compute pressure estimate for middle pressure from equation
            p_star = max(0, Ppvrs)
            if p_star <= pL :
                qL = 1.0	
            else:
                qL = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pL - 1))		

            SL = uL - aL*qL
            if p_star <= pR:
                qR = 1.0
            else:
                qR = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pR - 1))
            
            SR = uR + aR*qR

        elif(Wave_Estimate == 1): #Use Min-Mod Selection
            SL = min(uL-aL,uR-aR)
            SR = max(uL+aL,uR+aR)

        #Compute middle wave speed
        S_star = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR *(SR - uR))

        #DEBUG
        if(Debug_Flag == 1):
            print 'SL = ' + str(SL)
            print 'SR = ' + str(SR)
            print 'S* = ' + str(S_star)	

        if SL >= 0:
            #The flux is just the value of the left side flux
            FluxEstimates = CF.Compute_Cell_Fluxes(CellValues[0,:])
            return FluxEstimates

        elif SL <=0 and S_star >=0:
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


        elif S_star <= 0 and SR >=0:
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






class Domain(object):
    def __init__(self, input_parser):
        self.num_x = input_parser.user_input_data['Nx']
        self.x = None

    def initialize_coordinates(self):
        x_0 = input_parser.user_input_data['x0']
        x_L = input_parser.user_input_data['xL']
        
        dx = (x_L - x_0) / float(self.num_x)
        self.x = np.zeros(self.num_x)
        self.x[0] = dx / 2.0
        for i in range(1, self.num_x):
            self.x[i] = self.x[i-1] + dx

    def get_coordinates(self):
        return self.x

    def num_points(self):
        return len(self.x)


class BoundaryCondition(object):
    def __init__(self, u, p ,t, eos):
        self.eos = eos
        self.u = u
        self.p = p
        self.temperature = t

    def compute_conserved_values(self):
        bc_flux = np.zeros(3)
        bc_flux[0] = self.eos.density(self.p, self.temperature)
        bc_flux[1] = bc_flux[0] * self.u
        bc_flux[2] = bc_flux[0] * self.eos.energy(self.p, self.temperature, self.u)
        """
        print '\nBC Conserved for Left Boundary: \n'
        print 'rho = ' + str(bc_flux[0])
        print 'rho*u = ' + str(bc_flux[1])
        print 'rho*e = ' + str(bc_flux[2])
        """
        return bc_flux


class FlowInitializer(object):
    def __init__(self, input_parser):
        pass

    def initialize_solution(self, solution):
        raise NotImplementedError


class ShockTubeInitializer(FlowInitializer):
    def __init__(self, input_parser, eos):
        super(ShockTubeInitializer, self).__init__()
        self.eos = eos
        self.rho_1 = input_parser.user_input_data['rho1']
        self.u_1 = input_parser.user_input_data['u1']
        self.p_1 = input_parser.user_input_data['P1']
        
        self.rho_2 = .input_parser.user_input_data['rho2']
        self.u_2 = input_parser.user_input_data['u2']
        self.p_2 = input_parser.user_input_data['P2']
 
        self.shock_loc = input_parser.user_input_data['X_Loc_1']

    def initialize_solution(self, solution):
        print('Left and Right Initializations are:\n')
        print('Left Properties:\n')
        print('rho1 = ' + str(self.rho_1))
        print('rho1*u1 = ' + str(self.rho_1 * self.u_1))
        print('rho1*e1 = ' + str(self.rho_1 *(self.p_1/(self.rho_1*(gamma-1))) + 0.50*(u1*u1))

        print 'Right Properties:\n'
        print 'rho2 = ' + str(self.rho_2)
        print 'rho2*u2 = ' + str(self.rho_2 * self.u_2)
        print 'rho2*e2 = ' + str(self.rho_2 * (self.p_2 / (self.rho_2 * (gamma-1))) + 0.50*(u2*u2))

        print 'Shock is Initialized at:\n'
        print 'X_Shock = ' + str(shock_loc)

        #Loop over all cells and assign initial conditions
        for i in range(len(solution)):
            if self.solution.domain.x[i] <= self.shock_loc:
                solution[i]['rho'] = self.rho_1
                solution[i]['rho*u'] = self.rho_1 * self.u_1
                solution[i]['rho*e'] = self.rho_1 *self.eos.energy(self.p_1, self.rho_1) + 0.50*(self.u_1 * self.u_1)
            else:
                solution[i]['rho'] = self.rho_2
                solution[i]['rho*u'] = self.rho_2 * self.u_2
                solution[i]['rho*e'] = self.rho_2 * self.eos.energy(self.p_2, self.rho_2) + 0.50*(self.u_2 * self.u_2)



class Solution(object):
    import os
    def __init__(self, input_parser, domain, flow_initializer, solution_writer):
        self.input_parser = input_parser
        self.solution_writer = solution_writer
        self.time_integrator = time_integrator
        self.domain = domain
        self.flow_initializer = flow_initializer
        self.G = None
        self.initialize()

    def initialize(self):
        """
        # Creates 2D array for storing the conserved solution variables in the domain. 
        # Creates 2D array for storing the flux variables in the domain. 
        """
        # Array for holding the conserved cell variables(rho, rho*u, rho*E)
        self.G = [{'rho': 0, 'rho*u': 0, 'rho*e': 0} for i in range(self.domain.num_points)] #Solving 3 equations, so 3 columns
        
        self.flow_initialier.intialize_solution(self.G)

    def advance(self, dt):
           self.G = self.time_integrator(self.G, dt)

    def write_solution(self, iteration, current_time):
        """
        The purpose of this function is to output solution data to a file for post processing
        at a later time.
            Information: This function takes in the conserved variables array, G
        """
        # Before writing solution check to see if the output directory exists
        path = os.getcwd()
        if(os.path.isdir(path+'/output') == False): #Directory does not exist. Create it
            os.makedirs('output')

        OutPath = path+'/output'

        #Enter into output directory and print solution
        os.chdir(OutPath)
        
        #Create filename
        filename = 'solution_'+ str(iteration) +'_'+ str(current_time)

        #Open output file
        target = open(filename, 'w')

        #Compute primitives for output
        #rho, u, P, T
        prim = [self.compute_primitives(cell_value) for cell_value in self.G]

        #Output to file
        for i in range(self.domain.num_points()):
            outputString = str(x[i])+ '\t' +str(Prim[i,0])+'\t'+str(Prim[i,1])+'\t'+str(Prim[i,2])+'\t'+str(Prim[i,3])+'\n'
            target.write(outputString)

        target.close()


    def compute_primatives(self, cell_value):
        """
        The purpose of this function is to compute the values of the primitive variables from
        the flow field solution variables.
        """
        primitives = np.zeros((self.domain.num_points,4)) # rho, u, P, T    
        #Compute Density
        primitives[0] = cell_value['rho']

        #Compute velocity
        primitives[1] = cell_value['rho*u'] / cell_value['rho']

        #Compute Pressure
        primitives[2] = (gamma-1)*(CellValues[2] - 0.5*(CellValues[1]**2)/CellValues[0])

        #Compute Temperature: T = ( e - (1/2)u*u) * (gamma-1)/R
        primitives[3] = ( ( CellValues[2]/CellValues[0] ) - 0.5*(CellValues[1]/CellValues[0])**2)*(gamma-1)/R_gas

        return primatives    


class Solver(object):
    def __init__(self, initial_program_state):
        self.input_parser = initial_program_state['input_parser']
        self.flux_scheme = initial_program_state['flux_scheme']
        self.solution_writer = initial_program_state['solution_writer']
        self.time_integrator = initial_program_state['time_integrator']
        self.domain = initial_program_state['domain']

        self.num_timesteps = self.input_parser.user_input_data['Nt']
        self.dt = self.input_parser.user_input_data['dt']

    def initialize(self):
        self.domain.initialize_domain()


    def begin_timestepping(self):
        

        if(Real_Time_Plot == True):
            #Plotting Setup
            plt.ion()
            fig = plt.figure( 1,figsize = (10,10) )

        #Write Initial Domain Data
        self.solution.write_solution(0, 0*dt)
        for i in range(self.num_timesteps):
            self.solution.advance(self.dt)

            #Initialize printing variable
            Elapsed_Time = 0.0

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


            #Increment Elapsed Time
            Elapsed_Time = Elapsed_Time + dt

            if( Elapsed_Time >= Print_Time ):
                self.solution.write_solution(i+1, (i+1)*dt)
                Elapsed_Time = 0.0 #reset elapsed time

            print("Completed Timestep:\t%d\n"%(i+1))





