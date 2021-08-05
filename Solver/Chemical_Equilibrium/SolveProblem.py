"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Thur Oct 1 13:00:00 2020
----------------------------------------------------------------------
"""
from Solver.Chemical_Equilibrium.Equilibrate import equilibrate
from Solver.Chemical_Equilibrium.SolveProblemTP_TV import SolveProblemTP_TV



def SolveProblem(self, i):
    if any(self.PD.ProblemType.upper() == pt for pt in ['TP', 'TV']):
        self.PS.strP.append(SolveProblemTP_TV(self, self.PS.strR[i], self.PD.pP.Value, self.PD.TP.Value))
    else:
        if not i:
            self.PS.strP.append(equilibrate(self, self.PS.strR[i], self.PD.pP.Value))
        else:
            self.PS.strP.append(equilibrate(self, self.PS.strR[i], self.PD.pP.Value, self.PS.strP[i - 1]))    
    
    return self
    