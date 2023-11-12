# -*- coding: utf-8 -*-
"""
COMPUTE MORE STUFF NECESSARY TO INITIALIZE THE THERMOCHEMICAL CODE

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
"""
import numpy as np
from utils.databases.set_element_matrix import set_element_matrix

def Initialize_2(self):
    # Index gaseous and condensed species
    self = list_phase_species(self, self.S.LS)
    # Sort species: first gaseous species, secondly condensed species
    self.S = rearrange_species(self.S)
    # Stoichiometric Matrix
    self = Stoich_Matrix(self)
    return self


def list_phase_species(self, LS):
    self.S.ind_nswt = []
    self.S.ind_swt = []
    for ind, species in enumerate(LS):
        if not self.DB[species].swtCondensed:
            self.S.ind_nswt.append(ind)
        else:
            self.S.ind_swt.append(ind)
    # Number of gaseous species
    self.S.NG = len(self.S.ind_nswt)
    return self


def rearrange_species(self):
    self.LS = [self.LS[i] for i in self.ind_nswt + self.ind_swt]
    return self


def Stoich_Matrix(self):
    self.C.A0.Value = np.zeros((self.S.NS, self.E.NE))
    self.C.M0.Value = np.zeros((self.S.NS, 12))
    for i, species in enumerate(self.S.LS):
        txFormula = self.DB[species].txFormula
        self.DB[species].Element_matrix = set_element_matrix(
            txFormula, self.E.ElementsUpper)
        ind_Elements, atoms = (
            self.DB[species].Element_matrix[0, :], self.DB[species].Element_matrix[1, :])
        for ind_Element, atom in zip(ind_Elements, atoms):
            self.C.A0.Value[i, int(ind_Element)] = atom
        self.C.M0.Value[i, 9] = self.DB[species].swtCondensed

    self.C.N0.Value = self.C.M0.Value[:, [0, 9]]
    
    return self
