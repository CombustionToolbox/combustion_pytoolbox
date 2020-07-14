
# -*- coding: utf-8 -*-
"""
DISPLAY RESULTS
    
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Wen Jul 14 14:00:00 2020
----------------------------------------------------------------------
"""
import numpy as np
import pandas as pd

def displayResults(self, *args):
    ProblemType = self.PD.ProblemType
    mintol_display = self.C.mintol_display
    NameSpecies = self.S.NameSpecies
    
    strR = args[0]
    # strP = args[1]
    if len(args) == 3:
        str2 = args[2]
    print('-----------------------------------------------------------')
    print('PROBLEM TYPE: %s  |  EQUIVALENCE RATIO = %4.3f\n' % (ProblemType, strR.phi))
    
    props = ['T', 'p', 'r', 'h', 'e',
              's', 'cp', 'gamma']
    units = ['K', 'bar', 'kg/m3', 'kJ/kg', 'kJ/kg',
              'kJ/(kg-K)', 'kJ/(kg-K)', '-']
    hier_index = list(zip(props, units))
    hier_index = pd.MultiIndex.from_tuples(hier_index)
    data = np.array([
        [strR.T],
        [strR.p],
        [strR.rho],
        [strR.h/strR.mi],
        [strR.e/strR.mi],
        [strR.S],
        [strR.cP/strR.mi*1e-3],
        [strR.cP/strR.cV]
        ])
    df = pd.DataFrame(data, hier_index, ['REACTANTS'])
    df = df.round(3)
    df.index.names = ['PROPERTIES', 'UNITS']
    print(df)
    df = pd.DataFrame(strR.Xi, index=NameSpecies, columns=['Xi'])
    