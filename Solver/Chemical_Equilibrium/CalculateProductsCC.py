"""
COMPUTE CHEMICAL EQUILIBRIUM ASSUMING COMPLETE COMBUSTION
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Wen Jul 15 11:55:00 2020
----------------------------------------------------------------------
"""
import numpy as np
from Solver.Functions.SetSpecies import species_g0, equil_constant
from Solver.Chemical_Equilibrium.CalculatePhic import get_phic
from cmath import sqrt

def CalculateProductsCC(self, NatomE, phi, pP, TP):
    Elements, factor_c, Fuel, strThProp = [
        self.E.Elements, self.TN.factor_c, self.PD.Fuel, self.strThProp]
    R0 = self.C.R0  # [J/(K mol)]
    x, y, z = [NatomE[self.E.ind_C],
               NatomE[self.E.ind_H], NatomE[self.E.ind_O]]
    w, NHe, NAr = [NatomE[self.E.ind_N],
                   NatomE[self.E.ind_He], NatomE[self.E.ind_Ar]]  # Inerts

    NN2 = w/2
    NCgr = 0.

    FLAG_SOOT = False

    Ninerts = NHe + NAr
    phi_c = get_phic(self, Ninerts, TP, pP)
    # if Fuel.x and Fuel.x != Fuel.z:
    #     if Fuel.eps <= 0.5 and Fuel.eps > 1e-16:
    #         phi_c = (2 * (Fuel.x + Fuel.y/4 - Fuel.z/2)) / (-Fuel.z +
    #                 ((2 * Fuel.eps * Fuel.x + Fuel.x + Fuel.eps * Fuel.y/2)
    #                  / (1 + Fuel.eps)))  # C_x H_y O_z
    #     else:
    #         phi_c = 2/Fuel.x * (Fuel.x + Fuel.y/4 - Fuel.z/2)
    # else:
    #     phi_c = 1.1 * phi
    
    
    if phi <= 1.0: # Lean or stoichiometric mixtures
        NCO2 = x
        NCO = 0.
        NH2O = y/2
        NH2 = 0.
        NO2 =-x - y/4 + z/2
    else: # Rich mixtures
        NO2 = 0.
        if not x and y:
            NCO2 = 0.
            NCO = 0.
            NH2O = z
            NH2 = y/2 - z
        elif x and not y and phi < phi_c: # if there are only carbons (C) 
            NCO2 = -x + z
            NCO = 2*x - z
            NH2O = 0.
            NH2 = 0.
        elif phi < phi_c * factor_c:
            # General case of rich mixtures with hydrogens (H) and carbons (C)
            """         
            Equilibrium constant for the inverse wager-gas shift reaction
            
            CO2+H2 <-IV-> CO+H2O           
            """
            DG0 = (species_g0('CO', TP, strThProp)
                   + species_g0('H2O', TP, strThProp)
                   - species_g0 ('CO2', TP, strThProp)) * 1e3
            k4 = equil_constant(DG0, TP, R0)
            
            NCO = round((1/4) * (6*k4*x + k4*y - 2*k4*z - 4*x + 2*z
                            - np.sqrt(24*k4*x*z + 16*x**2 - 16*x*z - 16*k4*x**2 +
                                       4*k4**2*x*y - 8*k4**2 * x*z - 4 * k4**2 * y*z
                                       + 4*k4*y*z + 4*k4**2 * x**2 + k4**2 * y**2
                                       + 4*k4**2 * z**2 - 8*k4*z**2 +4*z**2)) / (k4-1), 14)
            NCO2 = x - NCO
            NH2O = -2*x + z + NCO
            NH2 = 2*x + y/2 - z - NCO
        elif phi >= phi_c * factor_c:
            # General case of rich mixtures with hydrogens (H) carbons (C) and soot
            """
            Equilibrium constant for the Boudouard reaction
            
            2CO <-VII-> CO2+C(gr)
            
            Equilibrium constant for the inverse wager-gas shift reaction
            
            CO2+H2 <-IV-> CO+H2O
            """
            DG0 = (species_g0('CO2', TP, strThProp) - 2*species_g0('CO', TP, strThProp)) * 1000
            k7 = equil_constant(DG0, TP, R0)
            DG0 = (species_g0('CO', TP, strThProp) + species_g0('H2O', TP, strThProp) - species_g0('CO2', TP, strThProp)) * 1000
            k4 = equil_constant(DG0, TP, R0)
            
            zeta = 1.0
            mu = k7 / zeta
            
            a0 = -2*k4 / zeta - k7*y + 2*k7*z
            a1 = 2*k7 + 4*k4*k7
            a2 = 4*k7**2 * zeta
            a3 = 2*k4*z / zeta
            
            NCO = np.real(-(a1/(3*a2))-(2**(1/3)*(-a1**2-3*a0*a2))/(3*a2*(-2*a1**3-9*a0*a1*a2+27*a2**2*a3+sqrt(-4*(a1**2+3*a0*a2)**3+(2*a1**3+9*a0*a1*a2-27*a2**2*a3)**2))**(1/3))+(-2*a1**3-9*a0*a1*a2+27*a2**2*a3+sqrt(-4*(a1**2+3*a0*a2)**3+(2*a1**3+9*a0*a1*a2-27*a2**2*a3)**2))**(1/3)/(3*2**(1/3)*a2))

            NCO2 = mu*NCO**2
            NCgr = x - NCO2 - NCO
            NH2O = z - 2*NCO2 - NCO
            NH2 = y/2 - NH2O
            
            FLAG_SOOT = True
    
    return [np.array([NCO2, NCO, NH2O, NH2, NO2, NN2, NHe, NAr, NCgr])[np.newaxis], phi_c, FLAG_SOOT]