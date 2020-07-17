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

def CalculateProductsIC(self, P, phi, pP, TP, vP, phi_c, FLAG_SOOT):
    E, S, PD, strThProp = [
        self.E, self.S, self.PD, self.TN, self.strThProp]
    M0, A0 = (self.C.M0.Value, self.C.A0.Value)
    R0TP = self.C.R0 * TP # [J/(mol)]
    
    it = 0
    itMax = 500
    t = True
    # Relaxation/iteration parameters
    relax = 0.00007385775 + (0.9854897 - 0.00007385775) / (1 + (TP/4058911)**1.817875)**658457.8
    # Number of moles of the major species in the product mixture under the
    # assumption of complete combustion (CC), denoted by subscript _0
    NCO2_0 = P_CC[S.idx_CO2, 0]
    NCO_0  = P_CC[S.idx_CO, 0]
    NH2O_0 = P_CC[S.idx_H2O, 0]
    NH2_0  = P_CC[S.idx_H2, 0]
    NO2_0  = P_CC[S.idx_O2, 0]
    NN2_0  = P_CC[S.idx_N2, 0]
    NCgr_0 = P_CC[S.idx_Cgr, 0]
    # Number of C, H, O, N, He, Ar-atoms in the product species
    NatomE = sum(P_CC[:, 0] * A0)
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
    NP = sum(P_CC[:, 0] * (1 - P_CC[:,9])) # Sum of num of moles of gases-(1-swt), with swt == condensed phase
    
    if 'P' in PD.ProblemType: # TP, HP, SP
        zeta = NP/pP
    elif 'V' in PD.ProblemType: # TV, EV, SV
        zeta = (vP * 1e-3) * 1e5 / R0TP 
    
    P_IC = M0.copy()    