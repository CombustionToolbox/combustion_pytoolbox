"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Thur Oct 1 13:00:00 2020
----------------------------------------------------------------------
"""
from Solver.Chemical_Equilibrium.Equilibrate import equilibrate
from Solver.Shocks_and_detonations.Shock_incident import shock_incident

def SolveProblem(self, i):
    if not any(self.PD.ProblemType.upper() == pt for pt in ['SHOCK_I', 'SHOCK_R', 'DET', 'DET_OVERDRIVEN']):
        if not i:
            self.PS.strP.append(equilibrate(self, self.PS.strR[i], self.PD.pP.Value))
        else:
            self.PS.strP.append(equilibrate(self, self.PS.strR[i], self.PD.pP.Value, self.PS.strP[i - 1]))
    elif self.PD.ProblemType.upper() == 'SHOCK_I':
        str1, str2 = shock_incident(self, self.PS.strR[i], self.PD.pR.Value, self.PD.TR.Value, self.PD.u1.Value)
        self.PS.strR.append(str1)
        self.PS.strP.append(str2)

    return self
    