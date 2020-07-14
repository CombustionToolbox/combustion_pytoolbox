# -*- coding: utf-8 -*-
"""
COMPUTE MORE STUFF NECESSARY TO INITIALIZE THE THERMOCHEMICAL CODE

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""


def Initialize_2(self):
    # Compute CHON equilibria for the minors products considered
    self = Compute_minors_species(self)
    # CH4 major species
    self.M = CH4_major(self.M)
    # OH major species
    self.M = OH_major(self.M)
    # Ask problem type
    self.PD.ProblemType = Ask_problem()
    return self


def Compute_minors_species(self):
    # First we eliminate from the minor species list those considered major
    # species in case the user has included any of them
    self.M.minors_products = [minor for minor in self.M.minors_products if not minor in self.S.List_fixed_Species]
    self.S.List_Compute_Species = self.S.List_fixed_Species + self.M.minors_products

    if self.Misc.FLAG_FIRST:
        self.M.L_minor = len(self.M.minors_products)
        if self.M.L_minor > 0:
            # Find index minors species
            self.M.ind_minor = [self.S.NameSpecies.index(minor) for minor in self.M.minors_products]
            # Properties of other minor species under consideration, which can
            # be written in the generic form C_alpha H_beta O_gamma N_omega
            self.C.alpha = [self.C.A0.Value[ind_minor, self.E.ind_C]
                            for ind_minor in self.M.ind_minor]
            self.C.beta = [self.C.A0.Value[ind_minor, self.E.ind_H]
                           for ind_minor in self.M.ind_minor]
            self.C.gamma = [self.C.A0.Value[ind_minor, self.E.ind_O]
                            for ind_minor in self.M.ind_minor]
            self.C.omega = [self.C.A0.Value[ind_minor, self.E.ind_N]
                            for ind_minor in self.M.ind_minor]

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
