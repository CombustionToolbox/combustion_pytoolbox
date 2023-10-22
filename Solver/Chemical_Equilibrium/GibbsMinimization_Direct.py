"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Tue Sep 15 09:37:00 2020
----------------------------------------------------------------------
"""

import numpy as np
import math
import pandas as pd # FOR CHECK
from numpy import log, exp
from scipy.optimize import fmin_slsqp
from Solver.Functions.SetSpecies import SetSpecies, species_g0


def equilibrium(self, pP, TP, strR):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0 = C.R0
    R0TP = R0 * TP # [J/(mol)]
    # Initialization
    NatomE = strR.NatomE
    NP_0 = 0.1
    NP = NP_0
    # N0 = N_CC
    N0[:, 0] = 0.1/S.N_Compute_Species
    # Gibbs free energy
    G0 = np.array([(species_g0(species, TP, strThProp)) / R0TP * 1e3 for species in S.List_Compute_Species])
    
    A0 = A0.transpose()
    def DGibbs0(nj):
        nonlocal G0
        NP = sum(nj)
        return np.sum(nj * (G0 + log(nj / NP)))
    def ec1(n):
        'equality constraint'
        return np.dot(A0, n) - NatomE.transpose()
    def ic1(n):
        '''inequality constraint
        all n>=0
        '''   
        return n    
    
    # Minimization of Gibbs free energy -> DG = 0
    N0[:, 0] = fmin_slsqp(DGibbs0, N0[:, 0], f_eqcons=ec1,f_ieqcons=ic1, iter=300, acc=C.tolN)
    
    DeltaNP = max(np.dot(A0, N0[:, 0]) - NatomE)
    return (N0, DeltaNP)