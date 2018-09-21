
class FaceFluxScheme(object):
    def __init__(self, eos, bc_left, bc_right):
        self.eos = eos
        self.bc_left = bc_left
        self.bc_right = bc_right

    def compute_face_fluxes(cell_values):
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
        cell_flux.append((density * velocity**2) + eos.pressure(cell_value))

        #Store flux for conservation of energy equation i.e. u*(rho*e + p)
        cell_flux.append(velocity * (density * energy + eos.pressure(cell_value)))

        return cell_flux

    
    def compute_cell_face_fluxes(G):
        """
        The purpose of this function is to compute the values of the cell face fluxes.

        Information: This function takes in a cell data Nx X 3 array that holds the 3 conserved
        variables for all cells in the domain.
            Output: A Nx X 3 X 2 array of the fluxes for mass, momentum,and energy at the west
            and east faces for each cell in the domain
        """


        #Holds to face fluxes on the left and right faces for each cell
        cell_face_fluxes = [{'left': 0, 'right': 0} for i in range(len(G))]
        
        #Holds the conserved variable data for cells to the left and right of a particular face
        adjacent_face_cells = {}

        #Compute and store left face flux for the first cell
        adjacent_face_cells['left_cell'] = self.bc_left.compute_bc_conserved()
        adjacent_face_cells['right_cell'] = G[0]

        self.compute_flux(adjacent_face_cells)
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
        pass

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



class ROEFluxScheme(FaceFluxScheme):
    def __init__(self, eos
        pass

    def compute_flux(cell_values):
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


