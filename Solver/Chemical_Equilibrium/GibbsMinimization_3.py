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
from Solver.Functions.SetSpecies import SetSpecies, species_g0

def remove_elements(NatomE, A0, tol):
    # Find zero sum elements
    ind_E0 = np.where(NatomE <= tol)
    ind_E0 = list(ind_E0[0])
    ind_A0_E0 = np.where(A0[:, ind_E0] > 0)
    return list(ind_A0_E0[0])

def temp_values(S, NatomE, tol):
    # List of indices with nonzero values
    temp_ind_E = np.where(NatomE > tol)
    temp_ind_E = list(temp_ind_E[0])
    temp_ind_nswt = S.ind_nswt[:] # copy by slicing
    temp_ind_swt = S.ind_swt[:] # copy by slicing
    temp_ind = temp_ind_nswt + temp_ind_swt
    temp_NE = len(temp_ind_E)
    temp_NS = S.NS 
    temp_ind_remove = []
    temp_LS = S.LS[:] # copy by slicing
    return (temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE, temp_NS, temp_LS, temp_ind_remove)

def remove_item(N0, zip1, zip2, ls1, ls2, cond0=True):
    ls0 = []
    for n, ind in zip(zip1, zip2):
        ls0.append(ind)
        if cond0:
            N0[ind, 0] = 0.
            if N0[ind, 1]:
                ls1.remove(ind)
            else:
                ls2.remove(ind)
                
    return (N0, ls0, ls1, ls2)

def update_temp(temp_ind_remove, temp_LS, temp_ind, temp_NS, LS0):
    if temp_ind_remove:
        for ind in temp_ind_remove:
            temp_LS.remove(LS0[ind])
        temp_ind = list(set(temp_ind) - set(temp_ind_remove))
        temp_NS = len(temp_ind)
        
    return (temp_LS, temp_ind, temp_NS)

def update_matrix(self):
    
    return self
def equilibrium(self, N_CC, phi, pP, TP, vP):
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0TP = C.R0 * TP # [J/(mol)]
    # Initialization
    NatomE = np.dot(N_CC[:, 0], A0)
    NP_0 = 0.1
    NP = NP_0
    
    it = 0
    itMax = 500
    # itMax = 50 + round(S.NS/2)
    SIZE = -log(C.tolN)
    e = 0.
    DeltaNP = 1.
    # Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
    # for the sum of elements whose value is <= tolN 
    ind_A0_E0 = remove_elements(NatomE, A0, self.C.tolN)
    # List of indices with nonzero values
    (temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E,
     temp_NE, temp_NS, temp_LS, temp_ind_remove) = temp_values(S, NatomE, self.C.tolN)
    # Remove species from the computed indeces list of gaseous and condensed species 
    # and append the indeces of species that we have to remove
    N0, temp_ind_remove, temp_ind_swt, temp_ind_nswt = remove_item(N0,N0[ind_A0_E0, 0], ind_A0_E0, temp_ind_swt, temp_ind_nswt)
    # Update temp values
    temp_LS, temp_ind, temp_NS = update_temp(temp_ind_remove, temp_LS, temp_ind, temp_NS, S.LS)
    # Initialize species vector N0 
    N0[temp_ind, 0] = 0.1/temp_NS
    # Dimensionless Standard Gibbs free energy 
    g0 = np.array([(species_g0(species, TP, strThProp)) * 1e3 for species in S.LS])
    G0RT = g0/R0TP
    # Construction of part of matrix A (complete)
    A11 = np.eye(temp_NS)
    A12 = -np.concatenate((A0[np.ix_(temp_ind, temp_ind_E)], np.ones(temp_NS).reshape(temp_NS, 1)), axis = 1)
    A1 = np.concatenate((A11, A12), axis=1)
    A22 = np.zeros((temp_NE + 1, temp_NE + 1))
    A0_T = A0.transpose()
            
    while DeltaNP > C.tolN and it < itMax:
        it += 1
        # Gibbs free energy
        G0RT[temp_ind_nswt] =  -(g0[temp_ind_nswt] / R0TP + log(N0[temp_ind_nswt, 0] / NP) + log(pP))
        # Construction of matrix A
        A21 = np.concatenate((N0[temp_ind, 0] * A0_T[np.ix_(temp_ind_E, temp_ind)], [N0[temp_ind, 0]]))
        A22[-1, -1] = -NP
        A = np.concatenate((A1, np.concatenate((A21, A22), axis=1)))
        # Construction of vector b            
        bi_0 = np.array([NatomE[E] - np.dot(N0[temp_ind, 0], A0[temp_ind, E]) for E in temp_ind_E])
        NP_0 = NP - sum(N0[temp_ind_nswt, 0])
        b = np.concatenate((G0RT[temp_ind], bi_0, np.array([NP_0])))
        # Solve of the linear system A*x = b
        x = np.linalg.solve(A, b)
        # Calculate correction factor
        e = []
        # sum_elements = np.dot(N0[temp_ind, 0], A0[np.ix_(temp_ind, temp_ind_E)])
        # BRATIO = min(sum_elements)/max(sum_elements)
        # if BRATIO < 1e-5:
        #     SIZE = log(1000)/BRATIO + log(1000) * 6.9077553
        # else:
        #     SIZE = -log(C.tolN) 
        for n, n_log_new in zip(N0[temp_ind, 0], x[0:temp_NS]):
            if log(n)/log(NP) <= -SIZE and n_log_new >= 0.:
                e.append(abs(-log(n/NP) - 9.2103404 / (n_log_new - x[-1])))
            else:
                e.append(min(2/max(5*abs(x[-1]), abs(n_log_new)), math.e**2))
        e = min(1, min(e))
           
        # Apply correction
        N0_log = log(N0[temp_ind, 0]) + e * x[0:temp_NS]
        NP_log = log(NP) + e * x[-1]
        # Apply antilog
        N0[temp_ind, :] = np.concatenate((exp(N0_log).reshape(len(temp_ind), 1), N0[temp_ind, 1].reshape(len(temp_ind), 1)), axis=1)
        
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
            temp_NS = len(temp_ind)
            A11 = np.eye(temp_NS)
            A12 = -np.concatenate((A0[np.ix_(temp_ind, temp_ind_E)], np.ones(temp_NS).reshape(temp_NS, 1)), axis = 1)
            A1 = np.concatenate((A11, A12), axis=1)
                    
        
        # print(f'\nit: {it}')
        # print(pd.DataFrame(N0[:, 0], index=np.array(S.LS)))
        NP = exp(NP_log) 
        DeltaN1 = max(np.array([n * abs(n_log) / NP for n, n_log in zip(N0[temp_ind, 0], x[0:temp_NS])]))
        DeltaN3 = NP_0 * abs(x[-1]) / NP
        DeltaNP = max(DeltaN1, DeltaN3) 
        # Deltab = [abs(bi - sum(N0[:, 0] * A0[:, i])) for i, bi in enumerate(x[S.NS:-1]) if bi > 1e-6]
        # print(Deltab)
    return (N0, DeltaNP)