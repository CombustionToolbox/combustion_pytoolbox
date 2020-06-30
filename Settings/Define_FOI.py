# -*- coding: utf-8 -*-
"""
DEFINE FUEL/OXIDIZER/INERT

Created on Tue Jun 30 16:28:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
from Solver.Functions.SetSpecies import SetSpecies
def Define_F(self, species, **moles):
    self.PD.S_Fuel = species
    self.PD.N_Fuel = [1.0]
    self.PD.R_Fuel = SetSpecies(self)
    return self

def Define_O(self, species, **moles):
    pass

def Define_I(self, species, **moles):
    pass