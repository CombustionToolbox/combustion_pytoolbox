"""
CALCULATE EQUILIBRIUM AT DEFINED T AND P (TP)
                        OR
CALCULATE EQUILIBRIUM AT DEFINED T AND CONSTANT V (TV)
INPUT:
    strR  = Prop. of reactives (phi,species,...)
    phi   = velocity upstream       [-]
    pP    = pressure of products    [bar]
    TP    = temperature of products [K]
OUTPUT:
    strP  = Prop. of products (phi,species,...)
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Wen Jul 15 11:55:00 2020
----------------------------------------------------------------------
"""
from .CalculateProductsCC import CalculateProductsCC
from Solver.Functions.SetSpecies import SetSpecies
from Solver.Functions.ComputeProperties import ComputeProperties

def SolveProblemTP_TV(self, strR, phi, pP, TP):
    moles_CC, phi_c, FLAG_SOOT = CalculateProductsCC(self, strR.NatomE, phi, TP)
    P = SetSpecies(self, self.S.List_Compute_Species, moles_CC.transpose(), TP)
    if self.PD.CompleteOrIncomplete.upper() == 'Incomplete':
        pass
    else:
        DeltaNP = 0.
    if self.PD.ProblemType[1] == 'P':
        strP = ComputeProperties(self, P, pP, TP)
    else:
        NP = sum(P[:, 0] * (1 - P[:, 9]))
        pP = (NP * TP * self.C.R0 / (strR.v/1e3)) / 1e5
        strP = ComputeProperties(self, P, pP, TP)
        
    strP.error_moles = DeltaNP
    
    return strP