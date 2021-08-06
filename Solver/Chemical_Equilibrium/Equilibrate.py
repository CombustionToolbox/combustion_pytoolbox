"""
CALCULATE CHEMICAL EQUILIBRIUM FOR A GIVEN TRANSFORMATION

INPUT:
    strR  = Prop. of reactives   (phi,species,...)
    pP    = pressure of products [bar]
OUTPUT:
    strP  = Prop. of products    (phi,species,...)
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Tue Aug 04 18:00:00 2021
----------------------------------------------------------------------
"""

import numpy as np
from Solver.Functions.Transformation import get_transformation
from Solver.Functions.SetSpecies import SetSpecies
from Solver.Functions.ComputeProperties import ComputeProperties
from Solver.Chemical_Equilibrium.GibbsMinimization import equilibrium

# from Solver.Chemical_Equilibrium.GibbsMinimization_numba import equilibrium
# from Solver.Chemical_Equilibrium.cython.GibbsMinimizationCython import equilibrium
# from Solver.Chemical_Equilibrium.GibbsMinimization_Soot_2 import equilibrium
# from Solver.Chemical_Equilibrium.GibbsMinimization_Reduced import equilibrium
# from Solver.Chemical_Equilibrium.GibbsMinimization_Direct import equilibrium # For checks

def equilibrate(self, strR, pP, strP=None):
    try:
        # get attribute of the specified transformations
        attr_name = get_attr_name(self)
        # compute initial guess
        guess = get_guess(self, strR, pP, attr_name, strP)
        # root finding: find the value x that satisfies f(x) = strP.xx(x) - strR.xx = 0
        x, ERR = root_finding(self, strR, pP, attr_name, guess)
        # compute properties
        strP = equilibrate_T(self, strR, pP, x)
        strP.error_problem = ERR
    except:
        print("An exception occurred: error Equilibrate.py")
        
    return strP

def equilibrate_T(self, strR, pP, TP):
    # Compute number of moles 
    N, DeltaNP = equilibrium(self, pP, TP, strR)
    # Compute properties of all species
    P = SetSpecies(self, self.S.LS, N[:, 0], TP)

    if self.PD.ProblemType[1] == 'P':
        strP = ComputeProperties(self, P, pP, TP)
    else:
        NP = sum(P[:, 0] * (1 - P[:, 9]))
        pP = (NP * TP * self.C.R0 / (strR.v/1e3)) / 1e5
        strP = ComputeProperties(self, P, pP, TP)
        
    strP.error_moles = DeltaNP
    
    return strP


def get_attr_name(self):
    if any(self.PD.ProblemType.upper() == pt for pt in ['TP', 'TV']):
        attr_name = 'T'
    elif self.PD.ProblemType.upper() == 'HP':
        attr_name = 'h'
    elif self.PD.ProblemType.upper() == 'EV':
        attr_name = 'e'
    elif any(self.PD.ProblemType.upper() == pt for pt in ['SP', 'SV']):
        attr_name = 'S'
        
    return attr_name


def get_guess(self, strR, pP, attr_name, strP):
    if any(self.PD.ProblemType.upper() == pt for pt in ['TP', 'TV']):
        return get_transformation(self, 'TP')
    elif strP:
        return strP.T
    return steff_guess(self, strR, pP, attr_name)

def get_point(self, x_vector, g_vector):
    return x_vector[1] - (x_vector[1] - x_vector[0]) / (g_vector[1] - g_vector[0]) * g_vector[1]


def get_point_aitken(x, g_vector):
    return x - (g_vector[0] - x)**2 / (g_vector[1] - 2*g_vector[0] + x)    


def get_gpoint(self, strR, pP, attr_name, x):
    strP = equilibrate_T(self, strR, pP, x)
    return (getattr(strP, attr_name) - getattr(strR, attr_name))


def get_partial_derivative(self, struct):
    ProblemType = self.PD.ProblemType.upper()
    if ProblemType == 'HP':
        value = struct.cP
    elif ProblemType == 'EV':
        value = struct.cV
    elif ProblemType == 'SP':
        value = struct.cP / struct.T
    elif ProblemType == 'SV':
        value = struct.cV / struct.T
    return value * 1e-3 # units


def get_ratio_newton(self, strR, pP, attr_name, x):
    strP = equilibrate_T(self, strR, pP, x)
    f = (getattr(strP, attr_name) - getattr(strR, attr_name))
    fprime = get_partial_derivative(self, strP)
    return (f, fprime)


def steff_guess(self, strR, pP, attr_name):
    x_l = strR.T + 50.
    x_r = 1500.
    
    #if np.isnan(get_gpoint(self, strR, pP, attr_name, x_l)):
    #    x_l = x_l + 100.
        
    g_l = get_gpoint(self, strR, pP, attr_name, x_l)
    g_r = get_gpoint(self, strR, pP, attr_name, x_r)
    
    if g_l * g_r > 0 and g_l < g_r:
        x = x_l - 50.
    elif g_l * g_r > 0 or (np.isnan(g_l) and np.isnan(g_r)):
        x = 2000.
    elif np.isnan(g_l) and not np.isnan(g_r):
        x = x_r - 100.
    elif not np.isnan(g_l) and np.isnan(g_r):
        x = x_l + 100.
    else:
        x = get_point(self, [x_l, x_r], [g_l, g_r])
        
    return x

            
def steff(self, strR, pP, attr_name, x, tol0=1e-3, itMax=30):
    """
        Steffenson method for finding roots
    """
    if any(self.PD.ProblemType.upper() == pt for pt in ['TP', 'TV']):
        return (get_transformation(self, 'TP'), 0)
    
    it = 0; g = 1.0; ERR = 1.0
    
    while (abs(ERR) > tol0 or abs(g) > tol0) and it < itMax:
        it = it + 1
        g = get_gpoint(self, strR, pP, attr_name, x)
        fx = abs(g - x)
        g_aux  = get_gpoint(self, strR, pP, attr_name, fx)
        fx2 = abs(g_aux - fx)
        if abs(fx2 - 2*fx + x) > tol0:
            x = get_point_aitken(x, [fx, fx2])
        else:
            x = fx
        ERR = abs(g_aux - g) / (1 + abs(g_aux))
    
    ERR = max(abs(ERR), abs(g))
    print_error(it, itMax, x, ERR)
    return (x, ERR)

def newton(self, strR, pP, attr_name, x0, tol0=1e-3, itMax=30):
    """
        Newton-Raphson method for finding roots
    """
    it = 0; ERR = 1.0
    
    while abs(ERR) > tol0 and it < itMax:
        it = it + 1
        f0, fprime0 = get_ratio_newton(self, strR, pP, attr_name, x0)
        x = abs(x0 - f0 / fprime0)
        
        f = get_gpoint(self, strR, pP, attr_name, x)
        ERR = max(abs((x - x0) / x), abs(f - f0))
        x0 = x
    
    
    print_error(it, itMax, x, ERR)
    return (x, ERR)


def root_finding(self, strR, pP, attr_name, x0, tol0=1e-3, itMax=30, method=steff):
    return method(self, strR, pP, attr_name, x0, tol0=1e-3, itMax=30)


def print_error(it, itMax, TP, ERR):
    if it > itMax:
        print('****************************\n')
        print('** Solution not converged **\n')
        print('** Temp  =  %4.2f         **\n', TP)
        print('** Error =  %4.2f%%       **\n', abs(ERR)*100)
        print('** It    =  %4.d          **\n', it)
        print('****************************\n')