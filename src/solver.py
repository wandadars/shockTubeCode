
class FaceFluxScheme(object):
    def __init__(self, eos):
        self.eos = eos

    def compute_cell_fluxes(cell_value):
        """
        The purpose of this function is to compute the value of the flux variables at the cell
        centers using the conserved variables for the cell.

        Information: This function takes in a cell data 1x3 array that holds the 3 conserved
        variables.
            Output: A 1x3 array of the fluxes for mass, momentum,and energy
        """
    

        #Allocate the array to hold the flux data
        cell_flux = np.zeros(3)

        #Store flux for conservation of mass equation i.e. rho*u
        cell_flux[0] = cell_value[1]

        #Store flux for conservation of momentum equation i.e. rho*u^2 + p
        density = cell_value[0]
        energy = cell_value[2] / cell_value[0]
        velocity = cell_value[1] / cell_value[0]
        cell_flux[1] = (cell_value[1]**2) / cell_value[0] + eos.pressure(density, energy, velocity) 

        #Store flux for conservation of energy equation i.e. u*(rho*e + p)
        cell_flux[2] = (cell_value[1] / cell_value[0]) * (cell_value[2] + eos.pressure(density, energy, velocity))

        return cell_flux

class HLLCFluxScheme(FaceFluxScheme):
    def __init__(self, eos



class Domain(object):
    def __init__(self, input_parser):
        self.input_parser = input_parser
        self.num_x = self.input_parser.user_input_data['Nx']
        self.x_0 = self.input_paresr.user_input_data['x0']
        self.x_L = self.input_paresr.user_input_data['xL']

    def initialize_coordinates(self):
        dx = (self.x_L - self.x_0) / float(self.num_x)
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
    def __init__(self, input_parser, solution):
        self.input_parser = input_parser

    def initialize_solution(self, solution):
        raise NotImplementedError


class ShockTubeInitializer(FlowInitializer):
    def __init__(self, input_parser, eos):
        super(ShockTubeInitializer, self).__init__(input_parser)
        self.eos = eos
        self.rho_1 = self.input_parser.user_input_data['rho1']
        self.u_1 = self.input_parser.user_input_data['u1']
        self.p_1 = self.input_parser.user_input_data['P1']
        
        self.rho_2 = self.input_parser.user_input_data['rho2']
        self.u_2 = self.input_parser.user_input_data['u2']
        self.p_2 = self.input_parser.user_input_data['P2']
 
        self.shock_loc = self.input_parser.user_input_data['X_Loc_1']

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
        for i in range(self.solution.domain.num_points()):
            if( self.solution.domain.x[i] <= X_Loc_1 ):
                solution[i,0] = self.rho_1
                solution[i,1] = self.rho_1 * self.u_1
                solution[i,2] = self.rho_1 *self.eos.energy(self.p_1, self.rho_1) + 0.50*(self.u_1 * self.u_1)
            else:
                solution[i,0] = self.rho_2
                solution[i,1] = self.rho_2 * self.u_2
                solution[i,2] = self.rho_2 * self.eos.energy(self.p_2, self.rho_2) + 0.50*(self.u_2 * self.u_2)



class Solution(object):
    import os
    def __init__(self, input_parser, domain, flow_initializer, solution_writer):
        self.input_parser = input_parser
        self.solution_writer = solution_writer
        self.domain = domain
        self.flow_initializer = flow_initializer
        self.G = None
        self.face_fluxes =None


    def initialize_solution(self):
        """
        # Creates 2D array for storing the conserved solution variables in the domain. 
        # Creates 2D array for storing the flux variables in the domain. 
        """
        # Array for holding the conserved cell variables(rho, rho*u, rho*E)
        self.G = np.zeros((self.domain.num_points, 3)) #Solving 3 equations, so 3 columns
        
        self.flow_initialier.intialize_solution(self.G)
        self.face_fluxes = np.zeros((self.domain.num_points, 3, 2))


    def compute_primatives(cell_value):
        """
        The purpose of this function is to compute the values of the primitive variables from
        the flow field solution variables.
        """

        primitives = np.zeros((self.domain.num_points,4)) # rho, u, P, T    

        #Compute Density
        primitives[0] = CellValues[0]

        #Compute velocity
        Primitives[1] = CellValues[1] / CellValues[0]

        #Compute Pressure
        Primitives[2] = (gamma-1)*(CellValues[2] - 0.5*(CellValues[1]**2)/CellValues[0])

        #Compute Temperature: T = ( e - (1/2)u*u) * (gamma-1)/R
        Primitives[3] = ( ( CellValues[:,2]/CellValues[:,0] ) - 0.5*(CellValues[:,1]/CellValues[:,0])**2)*(gamma-1)/R_gas

        return primatives    


    def write_ascii(self, iteration, current_time):
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
        Prim = np.zeros((Nx,4)) #rho, u,P,T
        Prim = compute_primitives(CellValues)

        #Output to file
        for i in range(self.domain.num_points()):
            outputString = str(x[i])+ '\t' +str(Prim[i,0])+'\t'+str(Prim[i,1])+'\t'+str(Prim[i,2])+'\t'+str(Prim[i,3])+'\n'
            target.write(outputString)

        target.close()




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
        for i in range(self.num_timesteps):
            self.solution.advance(self.dt)


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




