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
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Mon 26 9:30:00 2020
----------------------------------------------------------------------
"""

from .CalculateProductsCC import CalculateProductsCC
from .CalculateProductsIC import CalculateProductsIC
from Solver.Functions.SetSpecies import SetSpecies
from Solver.Functions.ComputeProperties import ComputeProperties
from Solver.Chemical_Equilibrium.SolveProblemTP_TV import SolveProblemTP_TV
# from Solver.Chemical_Equilibrium.GibbsMinimization_Soot_2 import equilibrium
# from Solver.Chemical_Equilibrium.GibbsMinimization_Reduced import equilibrium
# from Solver.Chemical_Equilibrium.GibbsMinimization_Direct import equilibrium # For checks

def SolveProblemHP(self, strR, phi, pP, TP):
    
    return strP