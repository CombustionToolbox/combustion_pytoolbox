"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Tue Sep 15 09:37:00 2020
----------------------------------------------------------------------
"""

import numpy as np
import pandas as pd # FOR CHECK
from numpy import log, exp
from Solver.Functions.SetSpecies import SetSpecies, species_g0

def equilibrium(self, N_CC, phi, pP, TP, vP):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0TP = C.R0 * TP # [J/(mol)]
    # Initialization
    NatomE = sum(N_CC[:, 0].reshape(S.N_Compute_Species, 1) * A0)
    NP_0 = sum(N_CC[:, 0] * (1.0 - N_CC[:,1])) # Sum of num of moles of gases-(1-swt), with swt == condensed phase
    NP_0 = 0.1
    NP = NP_0
    
    x0 = NatomE[E.ind_C]
    y0 = NatomE[E.ind_H]
    z0 = NatomE[E.ind_O]
    w0 = NatomE[E.ind_N]
    
    # N0 = N_CC
    N0[:, 0] = 0.1/S.N_Compute_Species
    it = 0
    itMax = 500
    SIZE = -log(10**-8)
    e = 0.
    DeltaNP = 1.
    # Gibbs free energy
    G0 = -np.array([(species_g0(species, TP, strThProp)) / R0TP * 1e3 for species in S.List_Compute_Species])
    while np.abs(DeltaNP / NP) > C.tolN and it < itMax:
        it += 1
        # Construction of matrix A
        A11 = np.eye(S.N_Compute_Species)
        A12 = -np.concatenate((A0, np.ones(S.N_Compute_Species).reshape(S.N_Compute_Species, 1)), axis = 1)
        A21 = np.concatenate((N0[:, 0] * A0.transpose(), [N0[:, 0]]))
        A22 = np.zeros((E.NE + 1, E.NE + 1))
        A22[-1, -1] = -NP
        A = np.concatenate((np.concatenate((A11, A12), axis=1), np.concatenate((A21, A22), axis=1)))
        # Construction of vector b
        bi_0 = np.array([NatomE[E] - sum(N0[:, 0] * A0[:, E]) for E in range(E.NE)])
        NP_0 = NP - sum(N0[:, 0] * (1.0 - N0[:, 1]))
        b = np.concatenate((G0, bi_0, np.array([NP_0])))
        # Solve of the linear system A*x = b
        x = np.linalg.solve(A, b)
        # Calculate correction factor
        e = []
        for n, n_log_new in zip(N0[:, 0], x[0:S.N_Compute_Species + 1]):
            if log(n)/log(NP) <= -SIZE and n_log_new >= 0.:
                e.append(abs(-log(n/NP) - 9.2103404 / (n_log_new - x[-1])))
            else:
                e.append(2/max(5*abs(x[-1]), abs(n_log_new)))
        e = min(1, min(e))
        DeltaNP = np.linalg.norm(np.array([NP - NP_0,
                                               x0 - sum(N0[:, 0] * A0[:, E.ind_C]),
                                               y0 - sum(N0[:, 0] * A0[:, E.ind_H]),
                                               z0 - sum(N0[:, 0] * A0[:, E.ind_O]),
                                               w0 - sum(N0[:, 0] * A0[:, E.ind_N])]))        
        # Apply correction
        N0_log = log(N0[:, 0]) + e * x[0:S.N_Compute_Species]
        NP_log = log(NP) + e * x[-1]
        # Apply antilog
        N0 = np.concatenate((exp(N0_log).reshape(S.N_Compute_Species, 1), N0[:, 1].reshape(S.N_Compute_Species, 1)), axis=1)
        # print(f'\nit: {it}')
        # print(pd.DataFrame(N0[:, 0], index=np.array(S.List_Compute_Species)))
        NP = exp(NP_log)
    return (N0, e)