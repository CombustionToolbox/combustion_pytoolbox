# -*- coding: utf-8 -*-
"""
Checks if the species is a reference element (e.g., 'C(gr)')

Created on Wen Jun 24 20:04:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""

def detect_location_of_phase_specifier(species):
    """
    Detect the location of the opening pharentesis of the phase identifier (if any) 
    """
    if '(' in species:
        n_open_parenthesis = species.rfind('(')
        n_close_parenthesis = species.rfind(')')
        if n_close_parenthesis == len(species)-1:
            n_open_parenthesis = n_open_parenthesis
        else:
            n_comma = species.find(',')
            if n_comma:
                if not (n_comma - n_close_parenthesis):
                    n_open_parenthesis = len(species)
            else:
                n_open_parenthesis = len(species)
    else:
        n_open_parenthesis = len(species)
    
    return n_open_parenthesis
                    