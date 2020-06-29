# -*- coding: utf-8 -*-
"""
COMPUTE MORE STUFF NECESSARY TO INITIALIZE THE THERMOCHEMICAL CODE

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
def Initialize_2(self):
    n_pass = []
    for i, minor in enumerate(self.M.minor_products):
        if not any(minor in self.S.List_fixed_Species):
            n_pass.append(i)
        