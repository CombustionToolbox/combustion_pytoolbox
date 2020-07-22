# -*- coding: utf-8 -*-
"""
COMPUTE MORE STUFF NECESSARY TO INITIALIZE THE THERMOCHEMICAL CODE

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np
from NASA_database.set_element_matrix import set_element_matrix

def Initialize_2(self):
    # List of species that we are going to compute
    self = Compute_Species(self)
    # Index fixed species
    self.S = Index_fixed_Species(self.S)
    # Stoichiometric Matrix
    self = Stoich_Matrix(self)
    # Compute CHON equilibria for the minors products considered
    self = Compute_minors_species(self)
    # CH4 major species
    self.M = CH4_major(self.M)
    # OH major species
    self.M = OH_major(self.M)
    # Ask problem type
    self.PD.ProblemType = Ask_problem()
    return self

def Compute_Species(self):
    # First we eliminate from the minor species list those considered major
    # species in case the user has included any of them
    self.M.minors_products = [minor for minor in self.M.minors_products if not minor in self.S.List_fixed_Species]
    self.S.List_Compute_Species = self.S.List_fixed_Species + self.M.minors_products
    self.S.N_Compute_Species = len(self.S.List_Compute_Species)
    return self

def Index_fixed_Species(self):
    self.ind_CO2 = self.List_Compute_Species.index('CO2')
    self.ind_CO = self.List_Compute_Species.index('CO')
    self.ind_H2O = self.List_Compute_Species.index('H2O')
    self.ind_H2 = self.List_Compute_Species.index('H2')
    self.ind_O2 = self.List_Compute_Species.index('O2')
    self.ind_N2 = self.List_Compute_Species.index('N2')
    self.ind_He = self.List_Compute_Species.index('He')
    self.ind_Ar = self.List_Compute_Species.index('Ar')
    self.ind_Cgr = self.List_Compute_Species.index('Cbgrb')

    self.ind_fixed = [self.ind_CO2, self.ind_CO, self.ind_H2O, self.ind_H2,
                      self.ind_O2, self.ind_N2, self.ind_He, self.ind_Ar, self.ind_Cgr]

    return self

def Stoich_Matrix(self):
    self.C.A0.Value = np.zeros((self.S.N_Compute_Species, self.E.NE))
    self.C.M0.Value = np.zeros((self.S.N_Compute_Species, 12))
    for i, species in enumerate(self.S.List_Compute_Species):
        txFormula = self.strThProp[species].txFormula
        self.strThProp[species].Element_matrix = set_element_matrix(
            txFormula, self.E.ElementsUpper)
        ind_Elements, atoms = (
            self.strThProp[species].Element_matrix[0, :], self.strThProp[species].Element_matrix[1, :])
        for ind_Element, atom in zip(ind_Elements, atoms):
            self.C.A0.Value[i, int(ind_Element)] = atom
        self.C.M0.Value[i, 9] = self.strThProp[species].swtCondensed

    self.C.N0.Value = self.C.M0.Value[:, [0, 9]]
    
    return self

def Compute_minors_species(self):
    if self.Misc.FLAG_FIRST:
        self.M.L_minor = len(self.M.minors_products)
        if self.M.L_minor > 0:
            # Find index minors species
            self.M.ind_minor = [self.S.List_Compute_Species.index(minor) for minor in self.M.minors_products]
            # Properties of other minor species under consideration, which can
            # be written in the generic form C_alpha H_beta O_gamma N_omega
            self.C.alpha = np.array([self.C.A0.Value[ind_minor, self.E.ind_C]
                            for ind_minor in self.M.ind_minor])
            self.C.beta = np.array([self.C.A0.Value[ind_minor, self.E.ind_H]
                           for ind_minor in self.M.ind_minor])
            self.C.gamma = np.array([self.C.A0.Value[ind_minor, self.E.ind_O]
                            for ind_minor in self.M.ind_minor])
            self.C.omega = np.array([self.C.A0.Value[ind_minor, self.E.ind_N]
                            for ind_minor in self.M.ind_minor])

            self.S.ind_all = self.S.ind_fixed + self.M.ind_minor
        else:
            self.S.ind_all = self.S.ind_fixed
    return self


def CH4_major(self):
    # Rich combustion of Hydrocarbons with oxygen produce CH4.
    # The nº of moles of the latter increase for richer mixtures
    # To estimate well the compounds of the final mixture we have to
    # recalculate the next species, which also play a pivotal role:
    # CH3, H, and CH.
    if 'CH4' in self.minors_products:
        self.major_CH4 = True
        self.ind_m_CH4 = self.minors_products.index('CH4')
        self.ind_m_CH3 = self.minors_products.index('CH3')
        self.ind_m_H = self.minors_products.index('H')
        self.ind_m_CH = self.minors_products.index('CH')
    else:
        self.major_CH4 = False
    return self


def OH_major(self):
    if 'OH' in self.minors_products:
        self.major_OH = True
        self.ind_m_OH = self.minors_products.index('OH')
    else:
        self.major_OH = False
    return self


def Ask_problem():
    ProblemType = 'TP'
    return ProblemType
