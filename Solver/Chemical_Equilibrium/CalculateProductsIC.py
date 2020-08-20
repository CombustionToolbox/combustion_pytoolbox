"""
COMPUTE CHEMICAL EQUILIBRIUM ASSUMING INCOMPLETE COMBUSTION
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Fri Jul 17 12:54:00 2020
----------------------------------------------------------------------
"""
import numpy as np
import pandas as pd
from Solver.Functions.SetSpecies import species_g0
import pandas as pd    
    
def print_Dataframe(N_IC, S, it):
    print(f'\nit: {it}')
    aux = np.array(S.List_Compute_Species)
    df = pd.DataFrame(N_IC[N_IC[:,0]>0.000001, 0], index=(aux[N_IC[:,0]>0.000001]))
    df.sort_values(0, ascending=False, inplace=True)
    print(df)   
     
def correction(x0, x1, TP):
    # Relaxation/iteration parameters
    relax = 0.00007385775 + (0.9854897 - 0.00007385775) / (1 + (TP/4058911)**1.817875)**658457.8
    return x0 + relax * (x1 - x0)

def indexation(N, N_old, n, index, TP):
    n_old = N_old[index, 0]
    if n_old:
        n = correction(n_old, n, TP)
    N[index, 0] = n
    return N

def correctionMajor(x0, x1, TP):
    if TP > 3000.0:
        x1 = correction(x0, x1, TP)
    if x1 < 0.0:
        x1 = 0.75 * x0
    return x1

def CalculateProductsIC(self, N_CC, phi, pP, TP, vP, phi_c, FLAG_SOOT):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (self.C.N0.Value, self.C.A0.Value)
    R0TP = self.C.R0 * TP # [J/(mol)]
    
    it = 0
    itMax = 500
    t = True
    
    # Number of moles of the major species in the product mixture under the
    # assumption of complete combustion (CC), denoted by subscript _0
    NCO2_0 = N_CC[S.ind_CO2, 0]
    NCO_0  = N_CC[S.ind_CO, 0]
    NH2O_0 = N_CC[S.ind_H2O, 0]
    NH2_0  = N_CC[S.ind_H2, 0]
    NO2_0  = N_CC[S.ind_O2, 0]
    NN2_0  = N_CC[S.ind_N2, 0]
    NCgr_0 = N_CC[S.ind_Cgr, 0]
    # Number of C, H, O, N, He, Ar-atoms in the product species
    NatomE = sum(N_CC[:, 0].reshape(
            N_CC[:, 0].size, 1) * A0)
    x = NatomE[E.ind_C]
    y = NatomE[E.ind_H]
    z = NatomE[E.ind_O]
    w = NatomE[E.ind_N]
    # Initial guess for the number of moles of the major species in the
    # product species under the assumption of incomplete combustion (IC)
    NCO2   = NCO2_0
    NCO    = NCO_0
    NH2O   = NH2O_0
    NH2    = NH2_0
    NO2    = NO2_0
    NN2    = NN2_0
    NCgr   = NCgr_0
    NHe    = NatomE[E.ind_He] 
    NAr    = NatomE[E.ind_Ar]
    # Initial guess for the overall number of moles of gaseous species in the
    # product mixture
    NP = sum(N_CC[:, 0] * (1.0 - N_CC[:,1])) # Sum of num of moles of gases-(1-swt), with swt == condensed phase
    
    if 'P' in PD.ProblemType: # TP, HP, SP
        zeta = NP/pP
    elif 'V' in PD.ProblemType: # TV, EV, SV
        zeta = (vP * 1e-3) * 1e5 / R0TP 
    # CALCULATION OF GIBBS FREE ENERGY, CONSTANTS EQUILIBRIUM AND OTHER
    g_CO2 = species_g0('CO2', TP, strThProp)
    g_CO = species_g0('CO', TP, strThProp)
    g_H2O = species_g0('H2O', TP, strThProp)
    DG0_I = (g_CO - g_CO2) * 1e3
    DG0_II = -g_H2O * 1e3
    k1 = np.exp(-DG0_I / R0TP)
    k2 = np.exp(-DG0_II / R0TP)
    
    if M.major_CH4:
        g_CH4 = species_g0('CH4', TP, strThProp)
        g_C2H2 = species_g0('C2H2_acetylene', TP, strThProp)
        g_C6H6 = species_g0('C6H6', TP, strThProp)
        DG0_VIII = (species_g0('CH3', TP, strThProp) - species_g0('H', TP, strThProp) - species_g0('CH4', TP, strThProp)) * 1e3
        
        k8 = np.exp(-DG0_VIII / R0TP)
        
        DG0_IX = -(g_C6H6 - 3*g_C2H2) * 1e3
        k9 = np.exp(-DG0_IX / R0TP)
        
        DG0_X = -(species_g0('H', TP, strThProp) - species_g0('CH', TP, strThProp)) * 1e3
        k10 = np.exp(-DG0_X / R0TP)
    
    if M.major_OH:
        g_OH = species_g0('OH', TP, strThProp)
        DG0_XI = g_OH * 1e3
        k11 = np.exp(-DG0_XI / R0TP)
    
    # Correction to the initial guess of oxygene moles and overall number of moles
    NO2 = NO2_0 + ((NH2O * k2 + NCO2 * k1) / 2)**(2/3) * zeta**(1/3)
    NP += NO2 - NO2_0
    if phi <= 1.0 and M.L_minor: # case lean-to-stoichiometric mixtures
        DNfactor_III = 1.0 - (C.beta + 2 * (C.gamma + C.omega)) / 4
        DG0_III = np.array([(species_g0(minor, TP, strThProp) - alpha * g_CO2
                    - (beta/2) * g_H2O) * 1e3 for minor, alpha, beta in zip(M.minors_products, C.alpha, C.beta)])
        k3 = np.exp(-DG0_III / R0TP)
    elif phi > 1.0: # case rich mixtures
        if not x and y and M.L_minor: # if there are only hydrogens (H)
            DG0_VI = np.array([(species_g0(minor, TP, strThProp) - alpha * g_CO2
                    - (gamma - 2*alpha) * g_H2O) * 1e3 for minor, alpha, gamma in zip(M.minors_products, C.alpha, C.gamma)])
            k6 = np.exp(-DG0_VI / R0TP)
            DNfactor_VI = 1.0 - C.alpha - (C.beta + C.omega) / 2
        # ..........
    # NESTED FUNCTIONS
    def incomplete_phi_1():
        nonlocal N_IC, it, NP, NCO2, NH2O, NH2, NO2, NN2 
        DeltaNP = 1.0
        while np.abs(DeltaNP / NP) > C.tolN and it < itMax:
            it += 1
            # Initialization of the product matrix for incomplete combustion (IC)
            NP_old = NP
            N_IC_old = N_IC
            N_IC = N0.copy()
            """
            In lean combustion the product mixture always contains O2, CO2,
            H2O, and N2 (in fuel-air combustion). The number of moles of
            CO and H2 can be calculated from the equilibrium condition for
            the reactions
            
                             CO2 <-I -> CO+(1/2) O2             [k1]
                             H2O <-II-> H2+(1/2) O2             [k2]
            
            For the remaining minor species, we calculate the number of
            moles from the equilibrium condition for the generalized
            reaction
            
            C.alpha * CO2 + (C.beta/2) * H2O + (C.gamma/2 - C.alpha - C.beta/4) * O2
            + (C.omega/2) * N2 <-III-> C_alpha * H_beta * O_gamma * N_omega [k3]
            
            Note that reactions 'i' and II are particular forms of the more
            general reaction III
            """
            NCO = NCO2 / NO2**(1/2) * zeta**(1/2) * k1
            N_IC = indexation(N_IC, N_IC_old, NCO, S.ind_CO, TP)
            NH2 = NH2O / NO2**(1/2) * zeta**(1/2) * k2
            N_IC = indexation(N_IC, N_IC_old, NH2, S.ind_H2, TP)
            
            # Determination of the number of moles of the minor species
            # from the equilibrium condition for the above reaction
            if M.L_minor:
                Ni = (k3 * NCO2**C.alpha * NH2O**(C.beta/2) * NN2**(C.omega/2)
                      * NO2**(C.gamma/2 - C.alpha - C.beta/4) * zeta**DNfactor_III)

                if M.major_OH:
                    Ni[M.ind_m_OH] = np.sqrt(NH2 * NO2 * k11 * zeta**(-3/2))
                Ni[Ni > NP_old] = 0.75 * Ni[Ni > NP_old]
                for ni, minor in zip(Ni, M.ind_minor):
                    N_IC = indexation(N_IC, N_IC_old, ni, minor, TP)
            
            # Correction of the number of moles of CO2, H2O, O2 and N2 from atom
            # conservation equations
            NCO2_old = NCO2
            NCO2 = NCO2_0 - sum(N_IC[:, 0] * A0[:, E.ind_C]) # C-atom conservation
            NCO2 = correctionMajor(NCO2_old, NCO2, TP)
            
            NH2O_old = NH2O
            NH2O = NH2O_0 - sum(N_IC[:, 0] * A0[:, E.ind_H]) / 2 # H-atom conservation
            NH2O = correctionMajor(NH2O_old, NH2O, TP)
            
            N_IC[[S.ind_CO2, S.ind_H2O], 0] = [NCO2, NH2O]
            
            NO2_old = NO2
            NO2 = NO2_0 + NCO2_0 + NCO_0/2 + NH2O_0/2 - sum(N_IC[:, 0] * A0[:, E.ind_O]) / 2 # O-atom conservation
            NO2 = correctionMajor(NO2_old, NO2, TP)
            
            NN2_old = NN2
            NN2 = NN2_0 - sum(N_IC[:, 0] * A0[:, E.ind_N]) / 2 # O-atom conservation
            NN2 = correctionMajor(NN2_old, NN2, TP)
            
            N_IC[[S.ind_O2, S.ind_N2], 0] = [NO2, NN2]
            N_IC[[S.ind_He, S.ind_Ar], 0] = [NHe, NAr]
            
            NP = sum(N_IC[:, 0] * (1.0 - N_IC[:, 1]))
            DeltaNP = np.linalg.norm(np.array([NP - NP_old,
                                               x - sum(N_IC[:, 0] * A0[:, E.ind_C]),
                                               y - sum(N_IC[:, 0] * A0[:, E.ind_H]),
                                               z - sum(N_IC[:, 0] * A0[:, E.ind_O]),
                                               w - sum(N_IC[:, 0] * A0[:, E.ind_N])]))
            
        return (N_IC, DeltaNP)    
    
    def incomplete_phi_2():
        nonlocal N_IC, it, NP, NCO2, NH2O, NH2, NO2, NN2 
        DeltaNP = 1.0
        while np.abs(DeltaNP / NP) > C.tolN and it < itMax:
            it += 1
            # Initialization of the product matrix for incomplete combustion (IC)
            NP_old = NP
            N_IC_old = N_IC
            N_IC = N0.copy()
            """
            In rich combustion without carbon atoms the product mixture
            contains H2O, H2 and N2 (in fuel-air combustion). The number
            of moles of these species must be calculated using the H, O
            and N-atom conservation equations
            
            For the minor species, we calculate the number of moles from
            the equilibrium condition for the generalized reaction
            
            C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
            
            with C.alpha = 0 if no carbon atoms are present.
            
            Determination of the number of moles of the minor species from
            the equilibrium condition for the above reaction
            """
            # Determination of the number of moles of the minor species
            if M.L_minor:
                Ni = (k6 * NCO2**C.alpha * NH2O**(C.gamma - 2*C.alpha) * NN2**(C.omega/2)
                      * NH2**(C.beta/2 - C.gamma + 2*C.alpha) * zeta**DNfactor_VI)

                if M.major_OH:
                    Ni[M.ind_m_OH] = np.sqrt(NH2 * NO2 * k11 * zeta**(-3/2))
                Ni[Ni > NP_old] = 0.75 * Ni[Ni > NP_old]
                for ni, minor in zip(Ni, M.ind_minor):
                    N_IC = indexation(N_IC, N_IC_old, ni, minor, TP)
            
            # Correction of the number of moles of O2
            NO2_old = NO2
            NO2 = zeta*(k2 * NH2O/NH2)**2
            if not NO2_old:
                NO2 = correction(NO2_old, NO2, TP)
            # Correction of the number of moles of H2O, H2 and N2 from
            # atom conservation equations
            NH2O_old = NH2O
            NH2O = NH2O_0 - 2*NO2 + sum(N_IC[:, 0] * A0[:, E.ind_O]) # O-atom conservation
            #NH2O = correctionMajor(NH2O_old, NH2O, TP)
            
            NH2_old = NH2
            NH2 = NH2_0 + 2*NO2 + sum(N_IC[:, 0] * A0[:, E.ind_O]) - sum(N_IC[:, 0] * A0[:, E.ind_H])/2 # H-atom conservation
            #NH2 = correctionMajor(NH2_old, NH2, TP)
            
            NN2_old = NN2
            NN2 = NN2_0 - sum(N_IC[:, 0] * A0[:, E.ind_N])/2 # O-atom conservation
            #NN2 = correctionMajor(NN2_old, NN2, TP)
            
            N_IC[[S.ind_H2O, S.ind_H2], 0] = [NH2O, NH2]
            N_IC[[S.ind_O2, S.ind_N2], 0] = [NO2, NN2]
            N_IC[[S.ind_He, S.ind_Ar], 0] = [NHe, NAr]
            
            NP = sum(N_IC[:, 0] * (1.0 - N_IC[:, 1]))
            DeltaNP = np.linalg.norm(np.array([NP - NP_old,
                                               y - sum(N_IC[:, 0] * A0[:, E.ind_H]),
                                               z - sum(N_IC[:, 0] * A0[:, E.ind_O]),
                                               w - sum(N_IC[:, 0] * A0[:, E.ind_N])]))
            
        return (N_IC, DeltaNP) 
            
    # CHEMICAL EQUILIBRIUM COMPUTATIONS
    N_IC = N0
    if phi <= 1.0: # case lean-to-stoichiometric mixtures
        return incomplete_phi_1()
    elif phi > 1.0 and not x and y: # case rich mixtures with only hydrogens (H)
        return incomplete_phi_2()
    elif (x and not y and phi < phi_c) and not FLAG_SOOT: # case rich mixtures with only carbons (C)
        return incomplete_phi_3()
    elif phi < phi_c * TN.factor_c and not FLAG_SOOT: # general case of rich mixtures with hydrogens (H) and carbons (C) and without soot
        return incomplete_phi_4()
    elif phi >= phi_c * TN.factor_c or FLAG_SOOT: # rich mixtures with soot
        if x and not y: # with only carbons (C)
            return incomplete_phi_5()
        return incomplete_phi_6() # general case
    
    ##### PRINT CONVERGENCE
    
    # return (N_IC, DeltaNP)
    
                
            