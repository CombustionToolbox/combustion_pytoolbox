"""
CALCULATION OF THE SOOT FORMATION EQUIVALENCE RATIO (PHI_C)
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Wen Sep 9 10:15:00 2020
----------------------------------------------------------------------
"""
import numpy as np
from Solver.Functions.SetSpecies import species_g0, equil_constant

def get_phic(self, Ninerts, TP, pP):
    x, y, z, w = [self.PD.Fuel.x, self.PD.Fuel.y, self.PD.Fuel.z, self.PD.Fuel.w]
    strThProp, proportion_N2_O2, R0 = [self.strThProp, self.PD.proportion_N2_O2, self.C.R0]
    if x and x != z:
        DG0 = (species_g0('CO2', TP, strThProp) - 2*species_g0('CO', TP, strThProp)) * 1e3
        k7 = equil_constant(DG0, TP, R0)
        DG0 = (species_g0('CO', TP, strThProp) + species_g0('H2O', TP, strThProp) - species_g0('CO2', TP, strThProp)) * 1e3
        k4 = equil_constant(DG0, TP, R0)
        
        phi_c0 = 2/(x - z) * (x + y/4 - z/2)
        phi_c = phi_c0
        it = 0; itMax = 20; tol = 1
        while tol > self.C.tolPhiSoot and it < itMax:
            it += 1
            NP = x + y/2 + Ninerts + w/2 + proportion_N2_O2/phi_c * (x + y/4 - z/2)
            zeta = NP/pP
            NCO = -((zeta - np.sqrt(zeta) * np.sqrt(4*k7*x + zeta))/(2*k7))
            NCO2 = x - NCO
            NH2O = y/(2*(1 + NCO/NCO2 * k4))
            
            phi_c_old = phi_c
            phi_c = 2*(x + y/4 - z/2) / (2*NCO2 + NH2O + NCO - z)
            tol = abs((phi_c - phi_c_old) / phi_c)
        
        return phi_c            
    return 1e5
    
    
    

