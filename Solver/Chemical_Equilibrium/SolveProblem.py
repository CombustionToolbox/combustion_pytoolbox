"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Thur Oct 1 13:00:00 2020
----------------------------------------------------------------------
"""

from Solver.Chemical_Equilibrium.SolveProblemTP_TV import SolveProblemTP_TV
from Solver.Chemical_Equilibrium.SolveProblemHP_EV import SolveProblemHP_EV, SolveProblemHP_EV_fast


def SolveProblem(self, i):
    if any(self.PD.ProblemType.upper() == pt for pt in ['TP', 'TV']):
        self.PS.strP.append(SolveProblemTP_TV(self, self.PS.strR[i], self.PD.pR.Value, self.PD.TP.Value))
    elif any(self.PD.ProblemType.upper() == pt for pt in ['HP', 'EV']):
        if not i:
            self.PS.strP.append(SolveProblemHP_EV(self, self.PS.strR[i], self.PD.pR.Value))
        else:
            self.PS.strP.append(SolveProblemHP_EV_fast(self, self.PS.strR[i], self.PD.pR.Value, self.PS.strP[i - 1]))
        
    
    return self
    