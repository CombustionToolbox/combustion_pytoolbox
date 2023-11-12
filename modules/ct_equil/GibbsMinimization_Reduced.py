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
from utils.SetSpecies import SetSpecies, species_g0

def equilibrium(self, N_CC, phi, pP, TP, vP):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0TP = C.R0 * TP # [J/(mol)]
    # Initialization
    NatomE = sum(N_CC[:, 0].reshape(S.NS, 1) * A0)
    NP_0 = sum(N_CC[:, 0] * (1.0 - N_CC[:,1])) # Sum of num of moles of gases-(1-swt), with swt == condensed phase
    NP_0 = 0.1
    NP = NP_0
    
    x0 = NatomE[E.ind_C]
    y0 = NatomE[E.ind_H]
    z0 = NatomE[E.ind_O]
    w0 = NatomE[E.ind_N]
    
    # N0 = N_CC
    N0[:, 0] = 0.1/S.NS
    it = 0
    itMax = 50 + round(S.NS/2)
    itMax = 300
    SIZE = -log(C.tolN)
    e = 0.
    DeltaNP = 1.
    # Dimensionless Standard Gibbs free energy 
    g0 = np.array([(species_g0(species, TP, strThProp)) * 1e3 for species in S.LS]) 
    G0RT = -g0 / R0TP
    # Construction of part of matrix A
    A11 = np.eye(S.NS)
    A12 = -np.concatenate((A0, np.ones(S.NS).reshape(S.NS, 1)), axis = 1)
    A1 = np.concatenate((A11, A12), axis=1)
    A22 = np.zeros((E.NE + 1, E.NE + 1))
    A0_T = A0.transpose()
    
    A11 = np.empty([E.NE, E.NE])
    while DeltaNP > 0.5*1e-5 and it < itMax:
        it += 1
        # Gibbs free energy
        G0RT[S.ind_nswt] =  -(g0[S.ind_nswt] / R0TP + log(N0[S.ind_nswt, 0] / NP) + log(pP))
        # Construction of matrix A
        for k in range(0, E.NE):
            for i in range(0, E.NE):
                A11[k, i] = sum(A0[:, k] * A0[:, i] * N0[:, 0])
        A12 = np.array([sum(A0[:, i] * N0[:, 0]) for i in range(E.NE)]).reshape(E.NE, 1)
        A1 = np.concatenate((A11, A12), axis=1)
        A2 = np.concatenate((A12.transpose(), [[sum(N0[:, 0] * (1.0 - N0[:, 1]) - NP)]]), axis=1)
        A = np.concatenate((A1, A2))
        
        # Construction of vector b
        bi_0 = np.array([NatomE[E] - sum(N0[:, 0] * A0[:, E] * (1.0 - G0RT)) for E in range(E.NE)])
        NP_0 = NP - sum(N0[:, 0] * (1.0 - N0[:, 1])) + sum(N0[:, 0] * G0RT)
        b = np.concatenate((bi_0, np.array([NP_0])))
        # Solve of the linear system A*x = b
        x = np.linalg.solve(A, b)
        # Compute log variation of moles
        ############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FINISH IT!
        # Relaxation parameter
        e = []
        # sum_elements = sum(N0[:, 0].reshape(S.NS, 1) * A0)
        # BRATIO = min(sum_elements)/max(sum_elements)
        # if BRATIO < 1e-5:
        #     SIZE = -log(C.tolN) # log(1000)/BRATIO + log(1000) * 6.9077553
        # else:
        #     SIZE = -log(C.tolN) 
        for n, n_log_new in zip(N0[:, 0], x[0:S.NS + 1]):
            # print(log(n)/log(NP), -SIZE, n_log_new)
            if log(n)/log(NP) <= -SIZE and n_log_new >= 0.:
                e.append(abs(-log(n/NP) - 9.2103404 / (n_log_new - x[-1])))
            else:
                # e.append(min(2/max(5*abs(x[-1]), abs(n_log_new)), math.e**2))
                e.append(2/max(5*abs(x[-1]), abs(n_log_new)))
        e = min(1, min(e))
           
        # Apply correction
        DeltaNi_log = np.array([x[-1] + sum(A0[j, :] * x[0:-1]) - G0RT[j] for j in range(S.NS)])
        N0_log = log(N0[:, 0]) + e * DeltaNi_log
        NP_log = log(NP) + e * x[-1]
        # Apply antilog
        N0 = np.concatenate((exp(N0_log).reshape(S.NS, 1), N0[:, 1].reshape(S.NS, 1)), axis=1)
        # for i, n in enumerate(N0[:, 0]):
        #     if log(n/NP) < -SIZE:
        #         N0[i, 0] = 0. 
        # print(f'\nit: {it}')
        # print(pd.DataFrame(N0[:, 0], index=np.array(S.LS)))
        NP = exp(NP_log)
        
        # DeltaNP = np.linalg.norm(np.array([x0 - sum(N0[:, 0] * A0[:, E.ind_C]),
        #                                    y0 - sum(N0[:, 0] * A0[:, E.ind_H]),
        #                                    z0 - sum(N0[:, 0] * A0[:, E.ind_O]),
        #                                    w0 - sum(N0[:, 0] * A0[:, E.ind_N])]))
        
        # DeltaNP = np.abs(DeltaNP / NP) 
        # DeltaN1 = max(np.array([n * abs(n_log) / NP * (1 - n_swt) for n, n_log, n_swt in zip(N0[:, 0], DeltaNi_log, N0[:, 1])]))
        # DeltaN2 = max(np.array([abs(n_log) / NP * n_swt for n, n_log, n_swt in zip(N0[:, 0], DeltaNi_log, N0[:, 1])]))
        # DeltaN3 = NP_0 * abs(x[-1]) / NP
        # DeltaNP = max(DeltaN1, DeltaN2, DeltaN3) 
        
        DeltaN1 = max(np.array([n * abs(n_log) / NP for n, n_log in zip(N0[:, 0], DeltaNi_log)]))
        DeltaN3 = NP_0 * abs(x[-1]) / NP
        DeltaNP = max(DeltaN1, DeltaN3)
        #Deltab = [abs(bi - sum(N0[:, 0] * A0[:, i])) for i, bi in enumerate(x[S.NS:-1]) if bi > 1e-6]
    return (N0, DeltaNP)