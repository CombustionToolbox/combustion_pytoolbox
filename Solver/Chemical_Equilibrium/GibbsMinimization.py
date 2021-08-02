"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Thur Oct 1 13:00:00 2020
----------------------------------------------------------------------
"""
import numpy as np
import math
import pandas as pd 
from numpy import log, exp
from memory_profiler import profile
from Solver.Functions.SetSpecies import SetSpecies, species_g0, get_tInterval

def remove_elements(NatomE, A0, tol):
    """ Find zero sum elements """
    ind_E0 = np.where(NatomE <= tol)
    ind_E0 = list(ind_E0[0])
    ind_A0_E0 = np.where(A0[:, ind_E0] > 0)
    return list(ind_A0_E0[0])

def temp_values(S, NatomE, tol):
    """ List of indices with nonzero values and lengths """
    temp_ind_E = np.where(NatomE > tol)
    temp_ind_E = list(temp_ind_E[0])
    temp_ind_nswt = S.ind_nswt[:] # copy by slicing
    temp_ind_swt = S.ind_swt[:] # copy by slicing
    temp_ind = temp_ind_nswt + temp_ind_swt
    temp_NE = len(temp_ind_E)
    temp_NS = S.NS 
    temp_ind_remove = []
    return (temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE, temp_NS, temp_ind_remove)

def remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE):
    """ Remove species from the computed indeces list of gaseous and condensed species 
        and append the indeces of species that we have to remove """
    ls0 = []
    for n, ind in zip(zip1, zip2):
        if log(n/NP) < -SIZE:
            ls0.append(ind)
            N0[ind, 0] = 0.
            if N0[ind, 1]:
                ls1.remove(ind)
            else:
                ls2.remove(ind)
                
    return (N0, ls0, ls1, ls2)

def update_temp(temp_ind, temp_NS, N0, zip1, zip2, ls1, ls2, NP, SIZE):
    """ Update temp items """
    N0, temp_ind_remove, temp_ind_swt, temp_ind_nswt = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE)
    if temp_ind_remove:
        temp_ind = list(set(temp_ind) - set(temp_ind_remove))
        temp_NS = len(temp_ind)
        
    return (temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS)

def update_matrix_A1(A0, A1, temp_NS, temp_NS0, temp_ind, temp_ind_E):
    """ Update stoichiometric submatrix A1 """
    if temp_NS < temp_NS0:
        A11 = np.eye(temp_NS)
        A12 = -np.concatenate((A0[np.ix_(temp_ind, temp_ind_E)], np.ones(temp_NS).reshape(temp_NS, 1)), axis = 1)
        return np.concatenate((A11, A12), axis=1), temp_NS0
    return A1, temp_NS0

def update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_E):
    """ Update stoichiometric submatrix A2 """
    A21 = np.concatenate((N0[temp_ind, 0] * A0_T[np.ix_(temp_ind_E, temp_ind)], [N0[temp_ind, 0]]))
    A22[-1, -1] = -NP
    return np.concatenate((A21, A22), axis=1)

def update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_E):
    """ Update stoichiometric matrix A """
    A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_E)
    return np.concatenate((A1, A2))

def update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, G0RT):
    """ Update coefficient vector b """
    bi_0 = np.array([NatomE[E] - np.dot(N0[temp_ind, 0], A0[temp_ind, E]) for E in temp_ind_E])
    NP_0 = NP - sum(N0[temp_ind_nswt, 0])
    return np.concatenate((-G0RT[temp_ind], bi_0, np.array([NP_0])))

def relax_factor(NP, zip1, zip2, DeltaNP, SIZE):
    """ Compute relaxation factor """
    e = []
    for n, n_log_new in zip(zip1, zip2):
            if log(n)/log(NP) <= -SIZE and n_log_new >= 0.:
                e.append(abs(-log(n/NP) - 9.2103404 / (n_log_new - DeltaNP)))
            else:
                e.append(min(2/max(5*abs(DeltaNP), abs(n_log_new)), math.e**2))
                
    return min(1, min(e))

def update_SIZE(N0, A0, temp_ind, temp_ind_E, tol):
    sum_elements = np.dot(N0[temp_ind, 0], A0[np.ix_(temp_ind, temp_ind_E)])
    BRATIO = min(sum_elements)/max(sum_elements)
    if BRATIO < 1e-5:
        return log(1000)/BRATIO + log(1000) * 6.9077553
    return -log(tol) 

def apply_antilog(N0, NP_log, temp_ind):
    N0[temp_ind, 0] = exp(N0[temp_ind, 0])
    NP = exp(NP_log)
    return (N0, NP)

def compute_STOP(NP_0, NP, DeltaNP, zip1, zip2):
    DeltaN1 = max(np.array([n * abs(n_log) / NP for n, n_log in zip(zip1, zip2)]))
    DeltaN3 = NP_0 * abs(DeltaNP) / NP
    # Deltab = [abs(bi - sum(N0[:, 0] * A0[:, i])) for i, bi in enumerate(x[S.NS:-1]) if bi > 1e-6]
    return max(DeltaN1, DeltaN3) 

def print_moles(N0, LS, it):
    """ Print number of moles of each species per iteration """
    print(f'\nit: {it}')
    print(pd.DataFrame(N0, index=np.array(LS)))

#@profile
def equilibrium(self, pP, TP, strR):
    """ Generalized Gibbs minimization method """
    E, S, C, M, PD, TN, strThProp = [self.E, self.S, self.C, self.M,
                                 self.PD, self.TN, self.strThProp]
    N0, A0 = (C.N0.Value, C.A0.Value)
    R0 = C.R0
    R0TP = R0 * TP # [J/(mol)]
    # Initialization
    NatomE = strR.NatomE
    NP_0 = 0.1
    NP = NP_0
    
    it = 0
    # itMax = 500
    itMax = 50 + round(S.NS/2)
    SIZE = -log(C.tolN)
    e = 0.
    STOP = 1.
    # Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
    # for the sum of elements whose value is <= tolN 
    ind_A0_E0 = remove_elements(NatomE, A0, self.C.tolN)
    # List of indices with nonzero values
    (temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E,
     temp_NE, temp_NS, temp_ind_remove) = temp_values(S, NatomE, self.C.tolN)
    # Update temp values
    temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS = update_temp(temp_ind, temp_NS, N0, N0[ind_A0_E0, 0], ind_A0_E0, temp_ind_swt, temp_ind_nswt, NP=self.C.tolN, SIZE=SIZE)
    temp_NS0 = temp_NS + 1;
    # Initialize species vector N0 
    N0[temp_ind, 0] = 0.1/temp_NS
    # Dimensionless Standard Gibbs free energy 
    g0 = np.array([(species_g0(species, TP, strThProp, get_tInterval(species, TP, self.strThProp), R0)) * 1e3 for species in S.LS])
    G0RT = g0/R0TP
    # Construction of part of matrix A (complete)
    A1, temp_NS0 = update_matrix_A1(A0, [], temp_NS, temp_NS0, temp_ind, temp_ind_E)
    A22 = np.zeros((temp_NE + 1, temp_NE + 1))
    A0_T = A0.transpose()
            
    while STOP > C.tolN and it < itMax:
        it += 1
        # Gibbs free energy
        G0RT[temp_ind_nswt] =  g0[temp_ind_nswt] / R0TP + log(N0[temp_ind_nswt, 0] / NP) + log(pP)
        # Construction of matrix A
        A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_E)
        # Construction of vector b            
        b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, G0RT)
        # Solve of the linear system A*x = b
        x = np.linalg.solve(A, b)
        # Calculate correction factor
        # update_SIZE(N0, A0, temp_ind, temp_ind_E, C.tolN)
        e = relax_factor(NP, N0[temp_ind, 0], x[0:temp_NS], x[-1], SIZE)   
        # Apply correction
        N0[temp_ind, 0] = log(N0[temp_ind, 0]) + e * x[0:temp_NS]
        NP_log = log(NP) + e * x[-1]
        # Apply antilog
        N0, NP = apply_antilog(N0, NP_log, temp_ind)
        # Update temp values in order to remove species with moles < tolerance
        temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS = update_temp(temp_ind, temp_NS, N0, N0[temp_ind, 0], temp_ind, temp_ind_swt, temp_ind_nswt, NP=NP, SIZE=SIZE)
        # Update matrix A
        A1, temp_NS0 = update_matrix_A1(A0, [], temp_NS, temp_NS0, temp_ind, temp_ind_E)
        # Print moles per iteration
        # print_moles(N0[:, 0], S.LS, it)
        # Compute STOP criteria
        STOP = compute_STOP(NP_0, NP, x[-1], N0[temp_ind, 0], x[0:temp_NS])
        
    return (N0, STOP)