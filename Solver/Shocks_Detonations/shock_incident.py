"""
CALCULATE INCIDENT SHOCK STATE AND POSTSTATE 

INPUT:
    str1  = Prop. of pre-state   (phi,species,...)
    pP    = pressure of products [bar]
OUTPUT:
    str2  = Prop. of post-state  (phi,species,...)
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Universidad Carlos III de Madrid
         
Last update Oct 13 2021
----------------------------------------------------------------------
"""

def Shock_incident(self, str1, u1, str2=None):
    try:
        # Abbreviations 
        C, TN = [self.C, self.TN]
        # Constants
        R0 = C.R0 # Universal gas constant [J/(mol-K)]
        # Assign values
        str1.u = u1       # velocity preshock [m/s] - laboratory fixed
        str1.w_shock = u1 # velocity preshock [m/s] - shock fixed
        # Miscelaneous
        it = 0
        itMax = TN.it_shocks
        STOP = 1.
        # Initial estimates of p2/p1 and T2/T1
        p2, T2, p2p1, T2T1 = get_guess(str1, str2, TN)
        # Loop
        while STOP > TN.tol_shocks and it < itMax:
            it += 1
            # Construction of the Jacobian matrix and vector b
            J, b = update_system(self, str1, p2, T2, R0)
            # Solve of the linear system J*x = b
            x = np.linalg.solve(J, b)
            # Calculate correction factor
            zeta = relax_factor(x)
            # Apply correction
            log_p2p1, log_T2T1 = apply_correction(x, p2p1, T2T1, zeta)
            # Apply antilog
            p2, T2 = apply_antilog(str1, log_p2p1, log_T2T1) # [Pa] and [K]
            # Update ratios
            p2p1 = p2 / (str1.p * 1e5)
            T2T1 = T2 / str1.T
            # Compute STOP criteria
            STOP = compute_STOP(x)

        # Check convergence
        print_convergence(STOP, TN.tol_shocks, T2)
        # Save state
        str2 = save_state(self, str1, T2, p2, STOP)

    except:
        print("An exception occurred: error Shock_incident.py")
        
    return str2


def get_guess(str1, str2, TN):
    if str2:
        V1 = 1/str1.rho # [m3/kg]
        V = V1/TN.volumeBoundRation

        p2p1 = 1 + (str1.rho * str1.u^2  / str1.p * (1 - V/V1)) * 1e-5
        T2T1 = p2p1 * V / V1

        p2 = p2p1 * str1.p * 1e5 # [Pa]
        T2 = T2T1 * str1.T       # [K]
           
    else:
        p2 = str2.p * 1e5 # [Pa]
        T2 = str2.T       # [K]

        p2p1 = p2 / (str1.p * 1e5)
        T2T1 = T2 / str1.T

    return (p2, T2, p2p1, T2T1)

def update_system(self, str1, p2, T2, R0):
    """ Update Jacobian matrix and vector b """
    r1 = str1.rho
    p1 = str1.p *1e5 # [Pa]
    T1 = str1.T
    u1 = str1.u
    W1 = str1.W * 1e-3 # [kg/mol]
    h1 = str1.h / str1.mi * 1e3 # [J/kg]
    # Calculate frozen state given T & p
    str2, r2, dVdT_p, dVdp_T = state(self, str1, T2, p2)
    
    W2 = str2.W * 1e-3
    h2 = str2.h / str2.mi * 1e3 # [J/kg]
    cP2 = str2.cP / str2.mi # [J/(K-kg)]
    
    alpha = (W1 * u1^2) / (R0 * T1)
    J1 = -r1/r2 * alpha * dVdp_T - p2 / p1
    J2 = -r1/r2 * alpha * dVdT_p
    b1 = p2/p1 - 1 + alpha * (r1/r2 - 1)
    
    J3 = -u1^2 / R0 * (r1/r2)^2 * dVdp_T + T2 / W2 * (dVdT_p - 1)
    J4 = -u1^2 / R0 * (r1/r2)^2 * dVdT_p - T2 * cP2 / R0
    b2 = (h2 - h1) / R0 - u1^2 / (2*R0) * (1 - (r1/r2)^2)
    
    J = [J1, J2, J3, J4]
    b = [b1, b2]

    return (J, b)

def state(self, str1, T, p):
    """ Calculate frozen state given T & p """
    self.PD.ProblemType = 'TP'
    p = p*1e-5 # [bar]
    str2 = equilibrate_T(self, str1, p, T)
    r2 = str2.rho
    dVdT_p = str2.dVdT_p
    dVdp_T = str2.dVdp_T

    return (str2, r2, dVdT_p, dVdp_T)

def relax_factor(x):
    """ Compute relaxation factor """
    factor = [0.40546511, 0.04879016]
    zeta = factor / abs(x)
    return min(1, min(zeta))

def apply_correction(x, p2p1, T2T1, zeta):
    """ Compute new estimates """
    log_p2p1 = log(p2p1) + zeta * x(1)
    log_T2T1 = log(T2T1) + zeta * x(2)

    return (log_p2p1, log_T2T1)

def apply_antilog(str1, log_p2p1, log_T2T1):
    """ compute p2 and T2 """
    p2 = np.exp(log_p2p1) * str1.p * 1e5 # [Pa]
    T2 = np.exp(log_T2T1) * str1.T

    return (p2, T2)

def compute_STOP(x):
    """ compute stop condition """
    return max(abs(x))


def save_state(self, str1, T2, p2, STOP):
    """ Save state """
    str2 = state(self, str1, T2, p2)
    str2.v_shock = str1.u * str1.rho / str2.rho
    str2.u = str1.u - str2.v_shock # velocity postshock [m/s] - laboratory fixed
    str2.error_problem = STOP

    return str2


def print_convergence(STOP, TOL, T):
    """ Print relative error """
    if STOP > TOL:
        print('***********************************************************\n')
        print('Convergence error: %.2f\n', STOP)

    if T > 2e4:
        print('***********************************************************\n')
        print('Validity of the next results compromise\n')
        print('Thermodynamic properties fitted to 20000 K\n')