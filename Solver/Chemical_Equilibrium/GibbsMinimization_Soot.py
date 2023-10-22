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
from Solver.Functions.SetSpecies import species_g0, get_tInterval

def equilibrium(self, pP, TP, strR):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0 = C.R0
    R0TP = R0 * TP # [J/(mol)]
    # Initialization
    NatomE = strR.NatomE
    NP_0 = sum(N0[S.ind_nswt, 0]) # Sum of number of moles of gases
    NP_0 = 0.1
    NP = NP_0
    
    # N0 = N_CC
    N0[:, 0] = 0.1/(S.NS-len(S.ind_swt))
    it = 0
    itMax = 500
    SIZE = -log(C.tolN)
    e = 0.
    DeltaNP = 1.
    # Dimensionless Standard Gibbs free energy 
    g0 = np.array([(species_g0(species, TP, strThProp, get_tInterval(species, TP, self.strThProp), R0)) * 1e3 for species in S.LS])
    G0RT = -g0 / R0TP
    # Construction of part of matrix A
    A11 = np.eye(S.NS)
    A11[S.ind_swt, S.ind_swt] = 0. # For condensed species
    A12 = -np.concatenate((A0, np.ones(S.NS).reshape(S.NS, 1)), axis = 1)
    A1 = np.concatenate((A11, A12), axis=1)
    A21 = np.concatenate((A0.transpose(), [np.zeros(S.NS)]))
    A22 = np.zeros((E.NE + 1, E.NE + 1))
    A0_T = A0.transpose()
    # List LS with nonzero values
    temp_ind_nswt = S.ind_nswt
    temp_ind_swt = S.ind_swt
    temp_ind = temp_ind_nswt + temp_ind_swt
    temp_NS = S.NS
    temp_NG = len(temp_ind_nswt)
    while DeltaNP > 0.5 * 1e-5 and it < itMax:
        it += 1
        # Gibbs free energy
        G0RT[temp_ind_nswt] =  -(g0[temp_ind_nswt] / R0TP + log(N0[temp_ind_nswt, 0] / NP) + log(pP))
        # Construction of matrix A
        A21[0:-1, temp_ind_nswt] = N0[temp_ind_nswt, 0] * A21[0:-1, temp_ind_nswt]
        A21[-1, temp_ind_nswt] = N0[temp_ind_nswt, 0]
        A22[-1, -1] = -NP
        A = np.concatenate((A1, np.concatenate((A21[:, temp_ind], A22), axis=1)))
        # Construction of vector b            
        bi_0 = np.array([NatomE[E] - np.dot(N0[temp_ind, 0], A0[temp_ind, E]) for E in range(E.NE)])
        NP_0 = NP - sum(N0[temp_ind_nswt, 0])
        b = np.concatenate((G0RT[temp_ind], bi_0, np.array([NP_0])))
        # Solve of the linear system A*x = b
        x = np.linalg.solve(A, b)
        # Calculate correction factor
        e = []
        sum_elements = np.dot(N0[temp_ind, 0], A0[temp_ind, :])
        BRATIO = min(sum_elements)/max(sum_elements)
        if BRATIO < 1e-5:
            SIZE = log(1000)/BRATIO + log(1000) * 6.9077553
        else:
            SIZE = -log(C.tolN) 
        for n, n_log_new in zip(N0[temp_ind, 0], x[temp_ind]):
            if log(n)/log(NP) <= -SIZE and n_log_new >= 0.:
                e.append(abs(-log(n/NP) - 9.2103404 / (n_log_new - x[-1])))
            else:
                e.append(min(2/max(5*abs(x[-1]), abs(n_log_new)), math.e**2))
        e = min(1, min(e))
           
        # Apply correction
        N0[temp_ind_nswt, 0] = log(N0[temp_ind_nswt, 0]) + e * x[0:temp_NG]
        N0[temp_ind_swt, 0] = N0[temp_ind_swt, 0] + e * x[temp_NG+1:temp_NS]
        NP = exp(log(NP) + e * x[-1])
        # Apply antilog
        N0[S.ind_nswt, 0] = exp(N0[S.ind_nswt, 0])
        temp_ind_remove = []
        for n, ind in zip(N0[temp_ind, 0], temp_ind):
            if log(n/NP) < -SIZE:
                N0[ind, 0] = 0.
                temp_ind_remove.append(ind)
                if N0[ind, 1]:
                    temp_ind_swt.remove(ind)
                else:
                    temp_ind_nswt.remove(ind)
        if temp_ind_remove:
            temp_ind = list(set(temp_ind) - set(temp_ind_remove))
            temp_NG = len(temp_ind_nswt)
            temp_NS = len(temp_ind)
            A11 = np.eye(temp_NS)
            A12 = -np.concatenate((A0[temp_ind, :], np.ones(temp_NS).reshape(temp_NS, 1)), axis = 1)
            A1 = np.concatenate((A11, A12), axis=1)
        # print(f'\nit: {it}')
        # print(pd.DataFrame(N0[:, 0], index=np.array(S.LS)))
        
        DeltaN1 = max(np.array([n * abs(n_log) / NP for n, n_log in zip(N0[temp_ind, 0], x[0:temp_NS])]))
        if S.ind_swt:
            DeltaN2 = max(np.array([abs(n_log) / NP for n, n_log in zip(N0[temp_ind_swt, 0], x[temp_ind_swt])]))
        else:
            DeltaN2 = 0.
        DeltaN3 = NP_0 * abs(x[-1]) / NP
        DeltaNP = max(DeltaN1, DeltaN2, DeltaN3) 
        
    return (N0, DeltaNP)