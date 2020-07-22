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
    NameSpecies = self.S.List_Compute_Species
    
    strR = args[0]
    strP = args[1]
    if len(args) == 3:
        str2 = args[2]
    print('-----------------------------------------------------------')
    print('PROBLEM TYPE: %s  |  EQUIVALENCE RATIO = %4.3f\n' % (ProblemType, strR.phi))
    # Print properties
    pd.options.display.float_format = '{:.3f}'.format
    props = ['T', 'p', 'r', 'h', 'e',
              's', 'cp', 'gamma']
    units = ['K', 'bar', 'kg/m3', 'kJ/kg', 'kJ/kg',
              'kJ/(kg-K)', 'kJ/(kg-K)', '-']
    hier_index = list(zip(props, units))
    hier_index = pd.MultiIndex.from_tuples(hier_index)
    data = np.array([
        [strR.T, strP.T],
        [strR.p, strP.p],
        [strR.rho, strP.rho],
        [strR.h/strR.mi, strP.h/strP.mi],
        [strR.e/strR.mi, strP.e/strP.mi],
        [strR.S, strP.S],
        [strR.cP/strR.mi *1e-3, strP.cP/strP.mi *1e-3] ,
        [strR.cP/strR.cV, strP.cP/strP.cV]
        ])
    df = pd.DataFrame(data, hier_index, ['REACTANTS', 'PRODUCTS'])
    df = df.round(3)
    df.index.names = ['PROPERTIES', 'UNITS']
    print(df, '\n')
    # Print reactants
    pd.options.display.float_format = '{:.4E}'.format
    pd.set_option('display.max_rows', None)
    
    df = pd.DataFrame(strR.Xi, index=NameSpecies, columns=['Xi [-]'])

    df.sort_values('Xi [-]', ascending=False, inplace=True)
    df.index.names = ['REACTANTS']
    print(df[df['Xi [-]']>0])
    
    df2 = pd.DataFrame([strR.Xi.sum()], index=['TOTAL    '], columns=[''])
    print(df2, '\n')
    # Print products
    df = pd.DataFrame(strP.Xi, index=NameSpecies, columns=['Xi [-]'])
    df.sort_values('Xi [-]', ascending=False, inplace=True)
    df.index.names = ['PRODUCTS ']
    print(df[df['Xi [-]']>mintol_display])
    
    df2 = pd.DataFrame([strP.Xi.sum()], index=['TOTAL    '], columns=[''])
    print(df2, '\n')
    
    