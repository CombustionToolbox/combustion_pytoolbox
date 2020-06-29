# -*- coding: utf-8 -*-
"""
Compute element matrix of species

Input:
    txFormula
Output:
    Element_matrix

Example for CO2

    txFormula = CO2
    Element_matrix = ( 5 7 )
                     ( 1 2 )

The species contains 1 atom of element 5 (C; 6 starting from 1) and
2 atoms of element 7 (O; 8 starting from 1) 

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid

Created on Wen Jun 24 20:04:00 2020
"""
import numpy as np

def set_element_matrix(txFormula, Elements):
    Element_matrix = np.array([]).reshape((2,0))
    for index in range(0, 5):
        Element = txFormula[index * 8: index * 8 + 2]
        if not '  ' in Element:
            Element = Element.replace(' ', '') # in case the element contains only one letter
            ind_Element = Elements.index(Element)
            num_Element = float(txFormula[index * 8 + 2: (index + 1) * 8])
            Element_matrix = np.hstack([Element_matrix, [[ind_Element], [num_Element]]])
    
    return Element_matrix
