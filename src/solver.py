import os


class Domain(object):
    def __init__(self, input_parser):
        self.num_x = input_parser.user_input_data['Nx']
        self.x = None
        self.initialize_coordinates()

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
        bc_values = {}
        bc_values['rho'] = self.eos.density_PT(self.p, self.temperature)
        bc_values['rho*u' = bc_values['rho'] * self.u
        bc_values['rho*e'] = bc_values['rho'] * self.eos.energy_PTU(self.p, self.temperature, self.u)
        """
        print '\nBC Conserved for Left Boundary: \n'
        print 'rho = ' + str(bc_flux[0])
        print 'rho*u = ' + str(bc_flux[1])
        print 'rho*e = ' + str(bc_flux[2])
        """
        return bc_values


class FlowInitializer(object):
    def __init__(self):
        pass

    def initialize_solution(self, solution):
        raise NotImplementedError


class ShockTubeInitializer(FlowInitializer):
    def __init__(self, input_parser, eos):
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
        print('rho1*e1 = ' + str(self.rho_1 *self.eos.energy_PRU(self.p_1, self.rho_1, self.u_1)) 

        print 'Right Properties:\n'
        print 'rho2 = ' + str(self.rho_2)
        print 'rho2*u2 = ' + str(self.rho_2 * self.u_2)
        print 'rho2*e2 = ' + str(self.rho_2 * self.eos.energy_PRU(self.p_2, self.rho_2, self.u_2))

        print 'Shock is Initialized at:\n'
        print 'X_Shock = ' + str(shock_loc)

        #Loop over all cells and assign initial conditions
        for i in range(len(solution.domain.x)):
            if solution.domain.x[i] <= self.shock_loc:
                solution.G[i]['rho'] = self.rho_1
                solution.G[i]['rho*u'] = self.rho_1 * self.u_1
                solution.G[i]['rho*e'] = self.rho_1 *self.eos.energy_PRU(self.p_1, self.rho_1, self.u_1)
            else:
                solution.G[i]['rho'] = self.rho_2
                solution.G[i]['rho*u'] = self.rho_2 * self.u_2
                solution.G[i]['rho*e'] = self.rho_2 * self.eos.energy_PRU(self.p_2, self.rho_2, self.u_2)


        print('\n Displaying Initialized Flowfield Vector\n')
        print('\t\t'+'X\t'+'rho'+'\t'+'rho*u'+'\t'+'rho*e\n')
        for i in range(len(solution.domain.x)):
            print('Cell ' + str(i + 1) + ':\t\t' + str(self.solution.domain.x[i]) + '\t'+ str(solution.G[i]['rho']) + '\t' + str(solution.G[i]['rho*u']) + '\t' + str(solution.G[i]['rho*e']) + '\n'



class Solution(object):
    def __init__(self, input_parser, flow_initializer, time_integrator, solution_writer):
        self.input_parser = input_parser
        self.solution_writer = solution_writer
        self.time_integrator = time_integrator
        self.domain = Domain(input_parser)
        self.flow_initializer = flow_initializer
        self.G = None
        self.initialize()

    def initialize(self):
        """
        Initializes the solution list based on whatever selection is made by the user
        """
        # list for holding the conserved cell variables(rho, rho*u, rho*E)
        self.G = [{'rho': 0, 'rho*u': 0, 'rho*e': 0} for i in range(self.domain.num_points)] #Solving 3 equations
        
        self.flow_initialier.intialize_solution(self.G)

    def advance(self, dt):
           self.G = self.time_integrator(self.G, dt)

    def write_solution(self, iteration, current_time):
        self.solution_writer.write_solution(self, iteration, current_time)


class SolutionWriter(object):
    def __init__(self, eos):
        self.eos = eos

    def write_solution(self):
        raise NotImplementedError

    def compute_primatives(self, cell_value):
        """
        The purpose of this function is to compute the values of the primitive variables from
        the flow field solution variables.
        """
        primitives = {}  # rho, u, P, T  

        #Compute Density
        primitives['rho'] = cell_value['rho']

        #Compute velocity
        primitives['u'] = cell_value['rho*u'] / cell_value['rho']

        #Compute Pressure
        primitives['P'] = self.eos.pressure(cell_value) 

        #Compute Temperature: T = ( e - (1/2)u*u) * (gamma-1)/R
        primitives['T'] = self.eos.temperature(cell_value) 

        return primatives    


class ASCIISolutionWriter(SolutionWriter):
    def __init__(self):

    def write_solution(self, solution, iteration, current_time):
        """
        The purpose of this function is to output solution data to a file for post processing
        at a later time.
        """
        # Before writing solution check to see if the output directory exists
        path = os.getcwd()
        if(os.path.isdir(path+'/output') == False): #Directory does not exist. Create it
            os.makedirs('output')

        OutPath = path+'/output'

        #Enter into output directory and print solution
        os.chdir(OutPath)
        
        #Create filename
        filename = 'solution_'+ str(iteration) + '_' + str(current_time)

        #Open output file
        target = open(filename, 'w')

        #Compute primitives for output
        #rho, u, P, T
        prim = [self.compute_primitives(cell_value) for cell_value in solution.G]

        #Output to file
        for i in range(solution.domain.num_points()):
            outputString = str(solution.domain.x[i]) + '\t' + str(Prim[i]['rho']) + '\t' + str(Prim[i]['u']) + '\t' + str(Prim[i]['P']) + '\t' + str(Prim[i]['T']) + '\n'
            target.write(outputString)

        target.close()



class Solver(object):
    def __init__(self, initial_program_state):
        self.input_parser = initial_program_state['input_parser']
        self.solution = initial_program_state['solution']

        self.num_timesteps = self.input_parser.user_input_data['Nt']
        self.dt = self.input_parser.user_input_data['dt']

    def begin_timestepping(self):
        if self.input_parser.user_input_data['real_time_plot'] == True:
            #Plotting Setup
            plt.ion()
            fig = plt.figure( 1,figsize = (10,10) )

        #Write Initial Domain Data
        self.solution.write_solution(0, 0*dt)
        for i in range(self.num_timesteps):
            self.solution.advance(self.dt)

            #Initialize printing variable
            Elapsed_Time = 0.0
            
            """
            if self.input_parser.user_input_data['real_time_plot'] == True:
                Prim = [compute_primitives(cell) for cell in self.solution.G]

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
                """

            #Increment Elapsed Time
            Elapsed_Time = Elapsed_Time + self.dt

            if Elapsed_Time >= float(self.input_parser.user_input_data['print_time']) :
                self.solution.write_solution(i + 1, (i + 1) *dt)
                Elapsed_Time = 0.0 #reset elapsed time

            print('Completed Timestep:\t%d\n'%(i+1))





