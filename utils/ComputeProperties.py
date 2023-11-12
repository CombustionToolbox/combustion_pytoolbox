# -*- coding: utf-8 -*-
"""
DEFINE FUEL/OXIDIZER/INERT

Created on Wen Jun 30 12:49:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
from utils.Mix import Mix

def ComputeProperties(self, SpeciesMatrix, p, T):
    mix = Mix(self, SpeciesMatrix, p, T)
    return mix
