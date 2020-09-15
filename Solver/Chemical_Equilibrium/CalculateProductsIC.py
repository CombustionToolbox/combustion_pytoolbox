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
import traceback
from cmath import sqrt
from Solver.Functions.SetSpecies import SetSpecies, species_g0
from .CalculateProductsCC import CalculateProductsCC

def print_Dataframe(N_IC, S, it):
    print(f'\nit: {it}')
    aux = np.array(S.List_Compute_Species)
    df = pd.DataFrame(N_IC[N_IC[:,0]>0.000001, 0], index=(aux[N_IC[:,0]>0.000001]))
    df.sort_values(0, ascending=False, inplace=True)
    print(df)   

def relax(TP):
    # Relaxation/iteration parameters
    return 0.00007385775 + (0.9854897 - 0.00007385775) / (1 + (TP/4058911)**1.817875)**658457.8

def correction(x0, x1, TP):
    return x0 + 0.1 * relax(TP) * (x1 - x0)

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

def getZeta(ProblemType, NP, pP, vP, R0TP):
    if 'P' in ProblemType: # TP, HP, SP
        return NP/pP
    elif 'V' in ProblemType: # TV, EV, SV
        return (vP * 1e-3) * 1e5 / R0TP
    
def CalculateProductsIC(self, N_CC, phi, pP, TP, vP, phi_c, FLAG_SOOT):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0TP = C.R0 * TP # [J/(mol)]
    
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
    
    zeta = getZeta(PD.ProblemType, NP, pP, vP, R0TP)
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
        
        DG0_X = -(species_g0('H', TP, strThProp) + species_g0('CH', TP, strThProp)) * 1e3
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
                               - (gamma - 2*alpha) * g_H2O) * 1e3 
                               for minor, alpha, gamma in zip(M.minors_products, C.alpha, C.gamma)])
            k6 = np.exp(-DG0_VI / R0TP)
            DNfactor_VI = 1.0 - C.alpha - (C.beta + C.omega)/2
        elif (x and not y and M.L_minor and phi < phi_c) and not FLAG_SOOT:    
            DG0_V = np.array([(species_g0(minor, TP, strThProp) - (gamma - alpha - beta/2) * g_CO2
                              - (beta/2) * g_H2O - (2*alpha - gamma + beta/2) * g_CO) * 1e3 
                              for minor, alpha, beta, gamma in zip(M.minors_products, C.alpha, C.beta, C.gamma)])
            k5 = np.exp(-DG0_V / R0TP)
            DNfactor_V = 1.0 - C.alpha - (C.beta + C.omega)/2
        elif phi < phi_c * TN.factor_c and not FLAG_SOOT: # general case of rich mixtures with hydrogens (H) and carbons (C)
            if M.L_minor:
                DG0_V = np.array([(species_g0(minor, TP, strThProp) - (gamma - alpha - beta/2) * g_CO2
                              - (beta/2) * g_H2O - (2*alpha - gamma + beta/2) * g_CO) * 1e3 
                              for minor, alpha, beta, gamma in zip(M.minors_products, C.alpha, C.beta, C.gamma)])
            
                k5 = np.exp(-DG0_V / R0TP)
                DNfactor_V = 1.0 - C.alpha - (C.beta + C.omega)/2
            DG0_IV = (g_CO + g_H2O - g_CO2) * 1e3
            k4 = np.exp(-DG0_IV / R0TP)
        elif phi >=  phi_c * TN.factor_c or FLAG_SOOT:
            if M.L_minor:
                DG0_V = np.array([(species_g0(minor, TP, strThProp) - (gamma - alpha - beta/2) * g_CO2
                              - (beta/2) * g_H2O - (2*alpha - gamma + beta/2) * g_CO) * 1e3 
                              for minor, alpha, beta, gamma in zip(M.minors_products, C.alpha, C.beta, C.gamma)])
            
                k5 = np.exp(-DG0_V / R0TP)
                DNfactor_V = 1.0 - C.alpha - (C.beta + C.omega)/2
            DG0_IV = (g_CO + g_H2O - g_CO2) * 1e3
            k4 = np.exp(-DG0_IV / R0TP)
            DG0_VII = (g_CO2 - 2*g_CO) * 1e3
            k7 = np.exp(-DG0_VII / R0TP)
            mu = k7/zeta
    # NESTED FUNCTIONS
    def incomplete_phi_1():
        nonlocal N_IC, it, NP, NCO2, NH2O, NO2, NN2 
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
            
        return (N_IC, DeltaNP, it)    
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
            NH2O = NH2O_0 - 2*NO2 - sum(N_IC[:, 0] * A0[:, E.ind_O]) # O-atom conservation
            NH2O = correctionMajor(NH2O_old, NH2O, TP)
            
            NH2_old = NH2
            NH2 = NH2_0 + 2*NO2 + sum(N_IC[:, 0] * A0[:, E.ind_O]) - sum(N_IC[:, 0] * A0[:, E.ind_H])/2 # H-atom conservation
            NH2 = correctionMajor(NH2_old, NH2, TP)
            
            NN2_old = NN2
            NN2 = NN2_0 - sum(N_IC[:, 0] * A0[:, E.ind_N])/2 # O-atom conservation
            NN2 = correctionMajor(NN2_old, NN2, TP)
            
            N_IC[[S.ind_H2O, S.ind_H2], 0] = [NH2O, NH2]
            N_IC[[S.ind_O2, S.ind_N2], 0] = [NO2, NN2]
            N_IC[[S.ind_He, S.ind_Ar], 0] = [NHe, NAr]
            
            NP = sum(N_IC[:, 0] * (1.0 - N_IC[:, 1]))
            DeltaNP = np.linalg.norm(np.array([NP - NP_old,
                                               y - sum(N_IC[:, 0] * A0[:, E.ind_H]),
                                               z - sum(N_IC[:, 0] * A0[:, E.ind_O]),
                                               w - sum(N_IC[:, 0] * A0[:, E.ind_N])]))
            
        return (N_IC, DeltaNP, it)
    def incomplete_phi_3():
        nonlocal N_IC, it, NP, NCO2, NCO, NH2O, NO2, NN2 
        DeltaNP = 1.0
        while np.abs(DeltaNP / NP) > C.tolN and it < itMax:
            it += 1
            # Initialization of the product matrix for incomplete combustion (IC)
            NP_old = NP
            N_IC_old = N_IC
            N_IC = N0.copy()
            """
            In rich combustion without hydrogen atoms the product mixture
            contains CO2, CO and N2 (in fuel-air combustion). The number
            of moles of these species must be calculated using the C, O
            and N-atom conservation equations
            
            For the minor species, we calculate the number of moles from
            the equilibrium condition for the generalized reaction
            
            (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
               +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
                            <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            
            with C.beta = 0 if no hydrogen atoms are present.
            
            Determination of the number of moles of the minor species from
            the equilibrium condition for the above reaction
            """
            # Determination of the number of moles of the minor species                    
            if M.L_minor:
                Ni = (k5 * NCO2**(C.gamma - C.alpha - C.beta/2) * NH2O**(C.beta/2) * NN2**(C.omega/2)
                      * NCO**(C.beta/2 - C.gamma + 2*C.alpha) * zeta**DNfactor_V)

                Ni[Ni > NP_old] = 0.75 * Ni[Ni > NP_old]
                for ni, minor in zip(Ni, M.ind_minor):
                    N_IC = indexation(N_IC, N_IC_old, ni, minor, TP)
            
            # Correction of the number of moles of O2
            NO2_old = NO2
            NO2 = zeta*(k1 * NCO2/NCO)**2
            if not NO2_old:
                NO2 = correction(NO2_old, NO2, TP)
            # Correction of the number of moles of CO2, CO and N2 from
            # atom conservation equations
            NCO2_old = NCO2
            NCO2 = NCO2_0 - 2*NO2 + 3*sum(N_IC[:, 0] * A0[:, E.ind_C]) - sum(N_IC[:, 0] * A0[:, E.ind_O]) # C-atom conservation
            NCO2 = correctionMajor(NCO2_old, NCO2, TP)
            
            NCO_old = NCO
            NCO = NCO_0 + 2*NO2 - 2*sum(N_IC[:, 0] * A0[:, E.ind_C]) + sum(N_IC[:, 0] * A0[:, E.ind_O]) # O-atom conservation
            NCO = correctionMajor(NCO_old, NCO, TP)
            
            NN2_old = NN2
            NN2 = NN2_0 - sum(N_IC[:, 0] * A0[:, E.ind_N])/2 # O-atom conservation
            NN2 = correctionMajor(NN2_old, NN2, TP)
            
            N_IC[[S.ind_CO2, S.ind_CO], 0] = [NCO2, NCO]
            N_IC[[S.ind_O2, S.ind_N2], 0] = [NO2, NN2]
            N_IC[[S.ind_He, S.ind_Ar], 0] = [NHe, NAr]
            
            NP = sum(N_IC[:, 0] * (1.0 - N_IC[:, 1]))
            DeltaNP = np.linalg.norm(np.array([NP - NP_old,
                                               x - sum(N_IC[:, 0] * A0[:, E.ind_C]),
                                               z - sum(N_IC[:, 0] * A0[:, E.ind_O]),
                                               w - sum(N_IC[:, 0] * A0[:, E.ind_N])]))
            
        return (N_IC, DeltaNP, it)
    def incomplete_phi_4():
        nonlocal N_IC, it, NP, NCO2, NCO, NH2O, NH2, NO2, NN2 
        DeltaNP = 1.0
        while np.abs(DeltaNP / NP) > C.tolN and it < itMax:
            it += 1
            # Initialization of the product matrix for incomplete combustion (IC)
            NP_old = NP
            N_IC_old = N_IC
            N_IC = N0.copy()
            """
            In rich combustion the product mixture always contains CO2, CO,
            H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            these species must be calculated combining the 4 atom
            conservation equations with the equilibrium condition for
            the inverse water-gas shift reaction
            
                           CO2+H2 <-IV-> H2O+CO             [k4]
            
            For the remaining minor species, we calculate the number of
            moles from the equilibrium condition for the generalized
            reactions
            
            (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
               +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
                           <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            
            or
            
            C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
        
            Reaction V is preferred for low hydrogen content (e.g., CO
            combustion), whereas reaction VI is preferred for low carbon
            content (e.g., H2 combustion). In general, we shall use
            reaction V, and leave reaction VI exclusively for H2
            combustion.
        
            To estimate the ammount of O2 in the product mixture we use the
            equilibrium condition for the reactions
        
            low Hydrogen:    CO2 <-I -> CO+(1/2) O2             [k1]
            low Carbon:      H2O <-II-> H2+(1/2) O2             [k2]
            
            Determination of the number of moles of the minor species from
            the equilibrium condition for the above reaction
            """
            # Determination of the number of moles of the minor species                    
            if M.L_minor:
                Ni = (k5 * NCO2**(C.gamma - C.alpha - C.beta/2) * NH2O**(C.beta/2) * NN2**(C.omega/2)
                      * NCO**(C.beta/2 - C.gamma + 2*C.alpha) * zeta**DNfactor_V)
                if M.major_CH4:
                    Ni[M.ind_m_CH4] = NH2 * Ni[M.ind_m_CH3] / (Ni[M.ind_m_H] * k8)
                Ni[Ni > NP_old] = 0.75 * Ni[Ni > NP_old]
                for ni, minor in zip(Ni, M.ind_minor):
                    N_IC = indexation(N_IC, N_IC_old, ni, minor, TP)
            # Check Ni
            # print('CHON: ', sum(N_IC[:,0] * A0[:,E.ind_C]), sum(N_IC[:,0] * A0[:,E.ind_H]), sum(N_IC[:,0] * A0[:,E.ind_O]), sum(N_IC[:,0] * A0[:,E.ind_N]))
            # Correction of the number of moles of O2
            NO2_old = NO2
            NO2 = zeta*(k1 * NCO2/NCO)**2
            if not NO2_old:
                NO2 = correction(NO2_old, NO2, TP)
            # Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            # atom conservation equations and equilibrium condition
            a = NCO2_0 + NCO_0 - sum(N_IC[:, 0] * A0[:, E.ind_C])
            b = NH2O_0 + NH2_0 - sum(N_IC[:, 0] * A0[:,E.ind_H])/2
            c = NCO_0 + NH2_0 + 2*NO2 - 2*sum(N_IC[:, 0] * A0[:,E.ind_C]) - sum(N_IC[:, 0] * A0[:,E.ind_H])/2 + sum(N_IC[:, 0] * A0[:,E.ind_O])
            d = b-c
            
            
            NCO2_old = NCO2
            NCO_old = NCO
            NH2O_old = NH2O
            NH2_old = NH2
            NN2_old = NN2
            
            NCO  = 0.5*((a + c) * k4 + d - np.sqrt((a - c)**2 * k4**2 + (2*(2*a*c + d*(a+c))) * k4 + d**2)) / (k4-1)
            NH2  = c - NCO
            NH2O = d + NCO
            NCO2 = a - NCO
            NN2  = NN2_0 - sum(N_IC[:,0] * A0[:,E.ind_N])/2 # N-atom conservation
            
            NCO2 = correctionMajor(NCO2_old, NCO2, TP)
            NCO = correctionMajor(NCO_old, NCO, TP)
            NH2O = correctionMajor(NH2O_old, NH2O, TP)
            NH2 = correctionMajor(NH2_old, NH2, TP)
            NN2 = correctionMajor(NN2_old, NN2, TP)
            
            N_IC[[S.ind_CO2, S.ind_CO], 0] = [NCO2, NCO]
            N_IC[[S.ind_H2O, S.ind_H2], 0] = [NH2O, NH2]
            N_IC[[S.ind_O2, S.ind_N2], 0] = [NO2, NN2]
            N_IC[[S.ind_He, S.ind_Ar], 0] = [NHe, NAr]
            
            NP = sum(N_IC[:, 0] * (1.0 - N_IC[:, 1]))
            DeltaNP = np.linalg.norm(np.array([NP - NP_old,
                                               x - sum(N_IC[:, 0] * A0[:, E.ind_C]),
                                               y - sum(N_IC[:, 0] * A0[:, E.ind_H]),
                                               z - sum(N_IC[:, 0] * A0[:, E.ind_O]),
                                               w - sum(N_IC[:, 0] * A0[:, E.ind_N])]))
            
        return (N_IC, DeltaNP, it)
    def incomplete_phi_5():
        try:
            pass
        except Exception as e:
            print('MODULE NOT IMPLEMENTED', repr(e))
    def incomplete_phi_6():
        nonlocal N_IC, it, NP, NCgr, NCO2, NCO, NH2O, NH2, NO2, NN2, zeta, NCgr_0, NCO2_0, NCO_0, NH2O_0, NH2_0, NN2_0, NO2_0 
        DeltaNP = 1.0
        while np.abs(DeltaNP / NP) > C.tolN and it < itMax:
            it += 1
            # Initialization of the product matrix for incomplete combustion (IC)
            NP_old = NP
            N_IC_old = N_IC
            N_IC = N0.copy()
            """
            In rich combustion the product mixture always contains CO2, CO,
            H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            these species must be calculated combining the 4 atom
            conservation equations with the equilibrium condition for
            the inverse water-gas shift reaction
        
                             CO2+H2 <-IV-> H2O+CO             [k4]
        
            For the remaining minor species, we calculate the number of
            moles from the equilibrium condition for the generalized
            reactions
            
            (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
               +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
                            <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            
            or
            
            C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
        
            Reaction V is preferred for low hydrogen content (e.g., CO
            combustion), whereas reaction VI is preferred for low carbon
            content (e.g., H2 combustion). In general, we shall use
            reaction V, and leave reaction VI exclusively for H2
            combustion.
            
            To estimate the ammount of O2 in the product mixture we use the
            equilibrium condition for the reactions
        
            low Hydrogen:    CO2 <-I -> CO+(1/2) O2             [k1]
            low Carbon:      H2O <-II-> H2+(1/2) O2             [k2]
            
            Determination of the number of moles of the minor species from
            the equilibrium condition for the above reaction
            """
            if t:
                # Determination of the number of moles of the minor species                    
                if M.L_minor:
                    Ni = (k5 * NCO2**(C.gamma - C.alpha - C.beta/2) * NH2O**(C.beta/2) * NN2**(C.omega/2)
                        * NCO**(C.beta/2 - C.gamma + 2*C.alpha) * zeta**DNfactor_V)
                    if M.major_CH4:
                        Ni[M.ind_m_CH4] = NH2 * Ni[M.ind_m_CH3] / (Ni[M.ind_m_H] * k8)
                    Ni[Ni > NP_old] = 0.75 * Ni[Ni > NP_old]
                    for ni, minor in zip(Ni, M.ind_minor):
                        N_IC = indexation(N_IC, N_IC_old, ni, minor, TP)
                # Check Ni
                # print('CHON: ', sum(N_IC[:,0] * A0[:,E.ind_C]), sum(N_IC[:,0] * A0[:,E.ind_H]), sum(N_IC[:,0] * A0[:,E.ind_O]), sum(N_IC[:,0] * A0[:,E.ind_N]))
                # Correction of the number of moles of O2
                NO2_old = NO2
                NO2 = zeta*(k1 * NCO2/NCO)**2
                if not NO2_old:
                    NO2 = correction(NO2_old, NO2, TP)
                # Correction of the number of moles of Cgr, CO, H2, CO2, H2O and N2 from
                # atom conservation equations and equilibrium condition
                a = NCO2_0 + NCO_0 + NCgr_0 - sum(N_IC[:, 0] * A0[:, E.ind_C])
                b = NH2O_0 + NH2_0 - sum(N_IC[:, 0] * A0[:, E.ind_H])/2
                c = 2*NCO2_0 + NCO_0 + NH2O_0 - 2*NO2 - sum(N_IC[:, 0] * A0[:, E.ind_O])
            
                NCO2_old = NCO2
                NCO_old = NCO
                NCgr_old = NCgr
                NH2O_old = NH2O
                NH2_old = NH2
                NN2_old = NN2
                
                NCgr = np.real(1/24*(24*a+(4*(2+k4))/(k4*mu)+(2*(1+1j*sqrt(3))*(-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu)))/(k4*(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3))+(2*1j*(1j+sqrt(3))*(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3))/(k4*mu**2)-(1/(6*k4**2*mu**3))*((2*(2+k4)*mu+((1+1j*sqrt(3))*mu**2*(-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu)))/(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3)+1j*(1j+sqrt(3))*(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3))**2))) 
                NCO2 = (1 + 2*a*mu - 2*mu*NCgr - sqrt(1 + 4*a*mu - 4*mu*NCgr))/(2*mu)
                
                if not np.real(NCO2):
                    NCgr = np.real(1/6*(6*a+(2+k4)/(k4*mu)+(4-2*k4+k4**2*(1-6*b*mu+6*c*mu))/(k4*(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3))+(1/(k4*mu**2))*((8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3))-(1/(6*k4**2*mu**3))*((2*mu+k4*mu-(mu**2*(-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu)))/(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3)+(8*mu**3-6*k4*mu**3-3*k4**2*mu**3+k4**3*mu**3-18*b*k4**2*mu**4-36*c*k4**2*mu**4-9*b*k4**3*mu**4+9*c*k4**3*mu**4+sqrt(mu**6*((-4+2*k4+k4**2*(-1+6*b*mu-6*c*mu))**3+(-8+6*k4+k4**3*(-1+9*b*mu-9*c*mu)+3*k4**2*(1+6*b*mu+12*c*mu))**2)))**(1/3))**2)))
                    NCO2 = (1 + 2*a*mu - 2*mu*NCgr - sqrt(1 + 4*a*mu - 4*mu*NCgr))/(2*mu)
                    
                NCO = (-1 + np.sqrt(1 + 4*a*mu - 4*mu*NCgr))/(2*mu)
                NH2O = c - 2*NCO2 - NCO
                NH2 = b - NH2O
                NN2  = NN2_0 - sum(N_IC[:,0] * A0[:,E.ind_N])/2 # N-atom conservation
                
                NCO2 = correctionMajor(NCO2_old, NCO2, TP)
                NCO = correctionMajor(NCO_old, NCO, TP)
                NH2O = correctionMajor(NH2O_old, NH2O, TP)
                NH2 = correctionMajor(NH2_old, NH2, TP)
                NN2 = correctionMajor(NN2_old, NN2, TP)
            
            if NCgr <= 1e-5 and it == 1:
                NCgr = 0.
                if t:
                    N_CC, phi_c, FLAG_SOOT = CalculateProductsCC(self, NatomE, phi, pP, TP)
                    P = SetSpecies(self, self.S.List_Compute_Species, N_CC[0, :], TP)
                    N_CC = P[:, [0, 9]]
                    
                    it = 0
                    
                    # Number of moles of the major species in the product mixture under the
                    # assumption of complete combustion (CC), denoted by subscript _0
                    NCO2_0 = N_CC[S.ind_CO2, 0]
                    NCO_0  = N_CC[S.ind_CO, 0]
                    NH2O_0 = N_CC[S.ind_H2O, 0]
                    NH2_0  = N_CC[S.ind_H2, 0]
                    NO2_0  = N_CC[S.ind_O2, 0]
                    NN2_0  = N_CC[S.ind_N2, 0]
                    NCgr_0 = N_CC[S.ind_Cgr, 0]
                    # Initial guess for the number of moles of the major species in the
                    # product species under the assumption of incomplete combustion (IC)
                    NCO2   = NCO2_0
                    NCO    = NCO_0
                    NH2O   = NH2O_0
                    NH2    = NH2_0
                    NO2    = NO2_0
                    NN2    = NN2_0
                    NCgr   = NCgr_0
                    # Initial guess for the overall number of moles of gaseous species in the
                    # product mixture
                    NP = sum(N_CC[:, 0] * (1.0 - N_CC[:,1])) # Sum of num of moles of gases-(1-swt), with swt == condensed phase
                    zeta = getZeta(PD.ProblemType, NP, pP, vP, R0TP)
                    # Correction to the initial guess of oxygene moles and overall number of moles
                    NO2 = NO2_0 + ((NH2O * k2 + NCO2 * k1) / 2)**(2/3) * zeta**(1/3)
                    NP += NO2 - NO2_0
                    

                # Determination of the number of moles of the minor species                    
                if M.L_minor:
                    Ni = (k5 * NCO2**(C.gamma - C.alpha - C.beta/2) * NH2O**(C.beta/2) * NN2**(C.omega/2)
                        * NCO**(C.beta/2 - C.gamma + 2*C.alpha) * zeta**DNfactor_V)
                    if M.major_CH4:
                        Ni[M.ind_m_CH4] = NH2 * Ni[M.ind_m_CH3] / (Ni[M.ind_m_H] * k8)
                    Ni[Ni > NP_old] = 0.75 * Ni[Ni > NP_old]
                    for ni, minor in zip(Ni, M.ind_minor):
                        N_IC = indexation(N_IC, N_IC_old, ni, minor, TP)
                # Check Ni
                # print('CHON: ', sum(N_IC[:,0] * A0[:,E.ind_C]), sum(N_IC[:,0] * A0[:,E.ind_H]), sum(N_IC[:,0] * A0[:,E.ind_O]), sum(N_IC[:,0] * A0[:,E.ind_N]))
                # Correction of the number of moles of O2
                NO2_old = NO2
                NO2 = zeta*(k1 * NCO2/NCO)**2
                if not NO2_old:
                    NO2 = correction(NO2_old, NO2, TP)
                # Correction of the number of moles of CO, H2, CO2, H2O and N2 from
                # atom conservation equations and equilibrium condition
                a = NCO2_0 + NCO_0 - sum(N_IC[:, 0] * A0[:, E.ind_C])
                b = NH2O_0 + NH2_0 - sum(N_IC[:, 0] * A0[:,E.ind_H])/2
                c = NCO_0 + NH2_0 + 2*(NO2 -NO2_0) - 2*sum(N_IC[:, 0] * A0[:,E.ind_C]) - sum(N_IC[:, 0] * A0[:,E.ind_H])/2 + sum(N_IC[:, 0] * A0[:,E.ind_O])
                d = b-c
                
                NCO2_old = NCO2
                NCO_old = NCO
                NH2O_old = NH2O
                NH2_old = NH2
                NN2_old = NN2
                
                NCO  = 0.5*((a + c) * k4 + d - np.sqrt((a - c)**2 * k4**2 + (2*(2*a*c + d*(a+c))) * k4 + d**2)) / (k4-1)
                NH2  = c - NCO
                NH2O = d + NCO
                NCO2 = a - NCO
                NN2  = NN2_0 - sum(N_IC[:,0] * A0[:,E.ind_N])/2 # N-atom conservation
                
                NCO2 = correctionMajor(NCO2_old, NCO2, TP)
                NH2O = correctionMajor(NH2O_old, NH2O, TP)
                NCO = correctionMajor(NCO_old, NCO, TP)
                NH2 = correctionMajor(NH2_old, NH2, TP)
                NN2 = correctionMajor(NN2_old, NN2, TP)
            
                    
            N_IC[[S.ind_CO2, S.ind_CO, S.ind_Cgr], 0] = [NCO2, NCO, NCgr]
            N_IC[[S.ind_H2O, S.ind_H2], 0] = [NH2O, NH2]
            N_IC[[S.ind_O2, S.ind_N2], 0] = [NO2, NN2]
            N_IC[[S.ind_He, S.ind_Ar], 0] = [NHe, NAr]
            
            NP = sum(N_IC[:, 0] * (1.0 - N_IC[:, 1]))
            DeltaNP = np.linalg.norm(np.array([NP - NP_old,
                                               x - sum(N_IC[:, 0] * A0[:, E.ind_C]),
                                               y - sum(N_IC[:, 0] * A0[:, E.ind_H]),
                                               z - sum(N_IC[:, 0] * A0[:, E.ind_O]),
                                               w - sum(N_IC[:, 0] * A0[:, E.ind_N])]))
            
        return (N_IC, DeltaNP, it)        
    # CHEMICAL EQUILIBRIUM COMPUTATIONS
    N_IC = N0
    if phi <= 1.0: # case lean-to-stoichiometric mixtures
        (N_IC, DeltaNP, it) = incomplete_phi_1()
        # print(f'Iterations {it}')
        return (N_IC, DeltaNP)
    elif phi > 1.0 and not x and y: # case rich mixtures with only hydrogens (H)
        (N_IC, DeltaNP, it) = incomplete_phi_2()
        # print(f'Iterations {it}')
        return (N_IC, DeltaNP)
    elif (x and not y and phi < phi_c) and not FLAG_SOOT: # case rich mixtures with only carbons (C)
        (N_IC, DeltaNP, it) = incomplete_phi_3()
        # print(f'Iterations {it}')
        return (N_IC, DeltaNP)
    elif phi < phi_c * TN.factor_c and not FLAG_SOOT: # general case of rich mixtures with hydrogens (H) and carbons (C) and without soot
        (N_IC, DeltaNP, it) = incomplete_phi_4()
        # print(f'Iterations {it}')
        return (N_IC, DeltaNP)
    elif phi >= phi_c * TN.factor_c or FLAG_SOOT: # rich mixtures with soot
        if x and not y: # with only carbons (C)
            return incomplete_phi_5()
        (N_IC, DeltaNP, it) = incomplete_phi_6()  # general case
        return (N_IC, DeltaNP)
    
    ##### PRINT CONVERGENCE
    
    # return (N_IC, DeltaNP)
    
                
            