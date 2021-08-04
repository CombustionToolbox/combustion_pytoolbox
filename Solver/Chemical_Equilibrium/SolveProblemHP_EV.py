"""
CALCULATE ADIABATIC T AND COMPOSITION AT CONSTANT P (HP)

INPUT:
    strR  = Prop. of reactives (phi,species,...)
    phi   = velocity upstream       [-]
    pP    = pressure of products    [bar]
    TP    = temperature of products [K]
OUTPUT:
    strP  = Prop. of products (phi,species,...)
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Tue Aug 03 11:00:00 2021
----------------------------------------------------------------------
"""

import numpy as np
from scipy.interpolate import interp1d
from Solver.Functions.SetSpecies import SetSpecies
from Solver.Functions.ComputeProperties import ComputeProperties
from Solver.Chemical_Equilibrium.SolveProblemTP_TV import SolveProblemTP_TV

def SolveProblemHP_EV(self, strR, pP):
    TP_l = 800.
    TP_r = 1500.
    attr_name = get_attr_name(self)
    strP = SolveProblemTP_TV(self, strR, pP, TP_l)
    if np.isnan(getattr(strP, attr_name)):
        TP_l = TP_l + 100.
        strP = SolveProblemTP_TV(self, strR, pP, TP_l)
    Q_l  = getattr(strP, attr_name) - getattr(strR, attr_name)
    strP = SolveProblemTP_TV(self, strR, pP, TP_r)
    Q_r  = getattr(strP, attr_name) - getattr(strR, attr_name)
    
    if Q_l * Q_r > 0 or (np.isnan(Q_l) and np.isnan(Q_r)):
        TP = 2500.
    elif abs(Q_l) < abs(Q_r) or abs(Q_l) >= abs(Q_r):
        TP = TP_r - (TP_r - TP_l) / (Q_r - Q_l) * Q_r
        strP = SolveProblemTP_TV(self, strR, pP, TP)
        Q  = getattr(strP, attr_name) - getattr(strR, attr_name)
        # TP = interp1d([Q_l, Q_r], [TP_l, TP_r])
    elif np.isnan(Q_l) or not np.isnan(Q_r):
        TP = TP_r - 100.
    elif not np.isnan(Q_l) or np.isnan(Q_r):
        TP = TP_l + 100.
    else:
        TP = TP_r - (TP_r - TP_l) / (Q_r - Q_l) * Q_r;
        strP = SolveProblemTP_TV(self, strR, pP, TP)
        Q  = getattr(strP, attr_name) - getattr(strR, attr_name)
        # TP = interp1d([Q_l, Q_r, Q], [TP_l, TP_r, TP])

    DeltaT = 1
    tol0 = 1e-10
    itMax = 30
    it = 0

    while (abs(DeltaT) > 1e-2 or abs(Q) > 1e-2) and it < itMax:
        it = it + 1
        strP = SolveProblemTP_TV(self, strR, pP, TP)
        Q  = getattr(strP, attr_name) - getattr(strR, attr_name)
        gx = abs(Q - TP)
        strP_aux = SolveProblemTP_TV(self, strR, pP, gx)
        Q_aux  = strP_aux.h - getattr(strR, attr_name)
        gx2 = abs(Q_aux - gx)
        if abs(gx2 - 2*gx + TP) > tol0:
            TP = TP - (gx - TP)**2 / (gx2 - 2*gx + TP)
        else:
            TP = gx;
        DeltaT = abs(Q_aux - Q) / (1 + abs(Q_aux))

    strP.error_problem = max(abs(DeltaT),abs(Q));

    return strP

def SolveProblemHP_EV_fast(self, strR, pP, strP):
    '''
    CALCULATE ADIABATIC T AND COMPOSITION AT CONSTANT P (HP)
                            OR
    CALCULATE EQUILIBRIUM COMPOSITION AT ADIABATIC T AND CONSTANT V (EV)

    INPUT:
        strR  = Prop. of reactives (phi,species,...)
        phi   = equivalence ratio       [-]
        pP    = pressure of products    [bar]
    OUTPUT:
        strP  = Prop. of products (phi,species,...)
    '''
    
    DeltaT = 1.
    tol0 = 1e-10
    itMax = 30
    it = 0
    TP = strP.T
    attr_name = get_attr_name(self)
    
    while (abs(DeltaT) > 1e-2 or abs(Q) > 1e-2) and it < itMax:
        it = it + 1
        strP = SolveProblemTP_TV(self, strR, pP, TP)
        Q  = getattr(strP, attr_name) - getattr(strR, attr_name)
        gx = abs(Q - TP)
        strP_aux = SolveProblemTP_TV(self, strR, pP, gx)
        Q_aux  = strP_aux.h - getattr(strR, attr_name)
        gx2 = abs(Q_aux - gx)
        if abs(gx2 - 2*gx + TP) > tol0:
            TP = TP - (gx - TP)**2 / (gx2 - 2*gx + TP)
        else:
            TP = gx
        DeltaT = abs(Q_aux - Q) / (1 + abs(Q_aux))
            
    return strP

def get_attr_name(self):
    if self.PD.ProblemType.upper() == 'HP':
        attr_name = 'h'
    else:
        attr_name = 'e'
        
    return attr_name