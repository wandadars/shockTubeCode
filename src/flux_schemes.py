
class FaceFluxScheme(object):
    def __init__(self, eos, bc_left, bc_right):
        self.eos = eos
        self.bc_left = bc_left
        self.bc_right = bc_right

    def compute_flux(cell_values):
        raise NotImplementedError

    def compute_states(self, adjacent_cell_values):
        raise NotImplementedError

    def print_cell_values(self, adjacent_cell_values):
        print 'Left Density = ' + str(adjacent_cell_values['left_cell']['rho'])
        print 'Left Momentum = ' + str(adjacent_cell_values['left_cell']['rho*u'])
        print 'Left Energy = ' + str(adjacent_cell_values['left_cell']['rho*e'])

        print 'Right Density = ' + str(adjacent_cell_values['right_cell']['rho'])
        print 'Right Momentum = ' + str(adjacent_cell_values['right_cell']['rho*u'])
        print 'Right Energy = ' + str(adjacent_cell_values['right_cell']['rho*e'])


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
        cell_flux.append((cell_value['rho'] * (cell_value['rho*u'] / cell_value['rho'])**2) + eos.pressure(cell_value))

        #Store flux for conservation of energy equation i.e. u*(rho*e + p)
        cell_flux.append((cell_value['rho*u'] / cell_value['rho']) * (cell_value['rho*e'] + eos.pressure(cell_value)))

        return cell_flux

    
    def compute_cell_face_fluxes(G):
        """
        The purpose of this function is to compute the values of the cell face fluxes.
        
        G: list of cells in domain
        """

        #Holds to face fluxes on the left and right faces for each cell
        cell_face_fluxes = [{'left': 0, 'right': 0} for i in range(len(G))]
        
        #Holds the conserved variable data for cells to the left and right of a particular face
        adjacent_face_cells = {}

        #Compute and store left face flux for the first cell
        adjacent_face_cells['left_cell'] = self.bc_left.compute_bc_conserved()
        adjacent_face_cells['right_cell'] = G[0]
        cell_face_fluxes[0]['left'] = self.compute_flux(adjacent_face_cells)

        #Compute and store right face flux for the first cell
        adjacent_face_cells['left_cell'] = G[0]
        adjacent_face_cells['right_cell'] = G[1]
        cell_face_fluxes[0]['right'] = self.compute_flux(adjacent_face_cells)
        
        #Inner cells have 2 face fluxes that must be determined, but the left
        #flux on an inner cell is the same as the right flux of the cell to the left
        #of the cell under consideration.
        for i in range(1, len(G) - 1):
            cell_face_fluxes[i]['left'] = cell_face_fluxes[i-1]['right']

            #Compute and store right face flux
            adjacent_face_cells['left_cell'] = G[i]
            adjacent_face_cells['left_cell'] = G[j+1]
            cell_face_fluxes[i]['right'] = self.compute_flux(adjacent_face_cells)



        #The Last Cell's left face flux is the same as the right flux from the cell to its left 
        cell_face_fluxes[-1]['left'] = cell_face_fluxes[-2]['right']
        
        #Compute and store right face flux due to boundary condition
        adjacent_face_cells['left_cell'] = G[-1]
        adjacent_face_cells['right_cell'] = self.bc_right.compute_bc_conserved()
        cell_face_fluxes[-1]['right'] = self.compute_flux(adjacent_face_cells)

        return cell_face_fluxes


class HLLCFluxScheme(FaceFluxScheme):
    def __init__(self, eos, bc_left, bc_right): 
        super(HLLCFluxScheme, self).__init__(eos, bc_left, bc_right)

    def compute_flux(cell_values):
        """
        # The purpose of this function is to compute an estimate for the flux across a cell
        # face using the HLLC flux approximation method.
        #
        # Information: This function takes in the states of the left and right cells and
        #              outputs the HLLC flux estimate for the face.
        """
        
        #cell_values is a dictionary with two elements. Each element maps to cell conserved variables
        #data to the left and right of a face.

        #Lowercase u variables correspond to x direction velocity. Uppercase U 
        #variables correspond to a 3 element vector containing conserved variables or
        #some quanity that is related to them.

        # Compute variables that are needed for the flux calculation and store in array
        states = self.compute_states(cell_states)
        
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
                qL = np.sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pL - 1))		

            SL = uL - aL * qL
            if p_star <= pR:
                qR = 1.0
            else:
                qR = sqrt( 1 + ((gamma+1)/(2*gamma))*(p_star/pR - 1))
            
            SR = uR + aR*qR

        elif Wave_Estimate == 1 : #Use Min-Mod Selection
            SL = min(uL-aL, uR-aR)
            SR = max(uL+aL, uR+aR)

        #Compute middle wave speed
        S_star = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / (rhoL * (SL - uL) - rhoR *(SR - uR))

        if Debug_Flag == 1 :
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


    def compute_states(self, adjacent_cell_values):
        #States is a dictionary that contains the relevant state information for the left
        #and right states that is computed from the conserved cell values i.e. primitive values.
        states = {} 

        #Store the velocities of the states 
        states['u_left'] = adjacent_cell_values['left_cell']['rho*u'] / adjacent_cell_values['left_cell']['rho'] #Left velocity
        states['u_right'] =  adjacent_cell_values['right_cell']['rho*u'] / adjacent_cell_values['right_cell']['rho'] #Right velocity
        
        #Compute speed of sound(SS)
        states['soundspeed_left'] = self.eos.soundspeed(adjacent_cell_values['left_cell']) #Left SS
        states['soundspeed_right'] = self.eos.soundspeed(adjacent_cell_values['right_cell']) #Right SS
        
        #Store density
        states['rho_left'] = adjacent_cell_values['left_cell']['rho']
        states['rho_right'] = adjacent_cell_values['right_cell']['rho']

        #Store Pressure
        states['pressure_left'] = self.eos.pressure(adjacent_cell_values['left_cell'])
        states['pressure_right'] = self.eos.pressure(adjacent_cell_values['right_cell'])

        #Store rho*e
        states['rho*e_left'] = adjacent_cell_values['left_cell']['rho*e']
        states['rho*e_right'] = adjacent_cell_values['right_cell']['rho*e']

        return states




class ROEFluxScheme(FaceFluxScheme):
    def __init__(self, eos, bc_left, bc_right): 
        super(ROEFluxScheme, self).__init__(eos, bc_left, bc_right)

    def compute_flux(self, adjacent_cell_values):
        # Compute variables that are needed for the flux calculation and store in array
        states = self.compute_states(adjacent_cell_values)

        #Store states into named variables to make calculations easier to read
        uL = states['u_left']
        uR = states['u_right']

        rhoL = states['rho_left']
        rhoR = states['rho_right']

        hL = states['enthalpy_left'] #Enthalpy.
        hR = states['enthalpy_right']

        #DEBUG
        if Debug_Flag == 1 :
            print '\nLeft & Right States Cell Face'
            print 'uL = ' + str(uL)
            print 'uR = ' + str(uR)
            print 'rhoL = ' + str(rhoL)
            print 'rhoR = ' + str(rhoR)
            print 'hL = ' + str(hL)
            print 'hR = ' + str(hR)


        #Compute density averaged quantities
        hBar = (np.sqrt(rhoL) * hL + np.sqrt(rhoR) * hR) / (np.sqrt(rhoL) + np.sqrt(rhoR))
        uBar = (np.sqrt(rhoL) * uL + np.sqrt(rhoR) * uR) / (np.sqrt(rhoL) + np.sqrt(rhoR))
        cBar = np.sqrt((gamma-1) * (hBar - 0.5*uBar*uBar))

        #Define characteristic speeds
        lambda1 = uBar - cBar
        lambda2 = uBar
        lambda3 = uBar + cBar

        #Create r_p vectors
        r_1 = []
        r_1.append(1)
        r_1.append(uBar - cBar)
        r_1.append(hBar - uBar*cBar)

        r_2 = []
        r_2.append(1)
        r_2.append(uBar)
        r_2.append(0.5*uBar*uBar)

        r_3 = []
        r_3.append(1)
        r_3.append(uBar + cBar)
        r_3.append(hBar + uBar*cBar)

        #Compute the r^p vectors for computing alpha_p
        r1 = []
        r1.append((uBar / (4*cBar)) * (2 + (gamma-1) * (uBar/cBar)))
        r1.append((-1.0 / (2.0*cBar)) * (1 + (gamma-1) * (uBar/cBar)))
        r1.append(0.5 * (gamma-1) * (1 / (cBar*cBar)))

        r2 = []
        r2.append(1 - 0.5 * (gamma-1) * ((uBar * uBar) / (cBar*cBar)))
        r2.append((gamma-1) * (uBar / (cBar*cBar)))
        r2.append(-(gamma-1) * (1.0 / (cBar*cBar)))

        r3 = []
        r3.append(-(uBar / (4*cBar)) * (2 - (gamma-1) * (uBar / cBar)))
        r3.append((1.0 / (2.0 * cBar)) * (1 - (gamma-1) * (uBar / cBar)))
        r3.append(0.5 * (gamma-1) * (1 / (cBar * cBar)))

        cell_diffs = [adjacent_cell_values[1][key] - adjacent_cell_values[0][key] for key in adjacent_cell_values[0].keys()]
        alpha1 = sum([r_val * diff for zip(r1, cell_diffs)]) 
        alpha2 = sum([r_val * diff for zip(r2, cell_diffs)]) 
        alpha3 = sum([r_val * diff for zip(r3, cell_diffs)]) 


        Term1 = abs(lambda1) * alpha1 * r_1
        Term2 = abs(lambda2) * alpha2 * r_2
        Term3 = abs(lambda3) * alpha3 * r_3
        FluxAverage = 0.5 * (self.compute_cell_flux(adjacent_cell_values[1]) + self.compute_cell_fluxes(adjacent_cell_values[0]))

        FluxEstimates = FluxAverage - 0.5 * (Term1 + Term2 + Term3)
        return FluxEstimates



    def compute_cell_states(self, adjacent_cell_values):
        """
        The purpose of this function is to compute a set of state variables for a cell to be
        used in the ROE flux calculation.
        """
        #States is a dictionary that contains the relevant state information for the left
        #and right states that is computed from the conserved cell values i.e. primitive values.
        states = {}

        #Store the velocities of the states
        states['u_left'] = adjacent_cell_values['left_cell']['rho*u'] / adjacent_cell_values['left_cell']['rho'] #Left velocity
        states['u_right'] = adjacent_cell_values['right_cell']['rho*u'] / adjacent_cell_values['right_cell']['rho'] #Right velocity

        #Store density
        states['rho_left'] = adjacent_cell_values['left_cell']['rho']
        states['rho_right'] = adjacent_cell_values['right_cell']['rho']

        #Store h(enthalpy)
        states['enthalpy_left'] = ( adjacent_cell_values['left_cell']['rho*e'] + self.eos.pressure(adjacent_cell_values['left_cell'])  ) / adjacent_cell_values['left_cell']['rho']
        states['enthalpy_right'] = ( adjacent_cell_values['right_cell']['rho*e'] + self.eos.pressure(adjacent_cell_values['right_cell'])  ) / adjacent_cell_values['right_cell']['rho']

        return states


