# -*- coding: utf-8 -*-
"""
DEFINE FUEL/OXIDIZER/INERT

Created on Tue Jun 30 16:28:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
from Solver.Functions.SetSpecies import SetSpecies
from Solver.Functions.ComputeProperties import ComputeProperties


def Define_F(self, species, **moles):
    self.PD.S_Fuel = species
    if species:  # not empty, i.e., there is/are fuel/s in the mixture
        self.PD.N_Fuel = [1.0]
        self.PD.R_Fuel = SetSpecies(self)
        self.PS.strR_Fuel = ComputeProperties(
            self, self.PD.R_Fuel, self.PD.pR.Value, self.PD.TR.Value)
        self.PD.Fuel.x = self.PS.strR_Fuel.NatomE[self.E.ind_C]
        self.PD.Fuel.y = self.PS.strR_Fuel.NatomE[self.E.ind_H]
        self.PD.Fuel.z = self.PS.strR_Fuel.NatomE[self.E.ind_O]
        self.PD.phi.t = self.PD.Fuel.x + self.PD.Fuel.y / 4 - self.PD.Fuel.z / 2
    return self


def Define_O(self, species, i, **moles):
    self.PD.S_Oxidizer = species
    if species:  # not empty, i.e., there is/are oxidizer/s in the mixture
        if moles:
            self.PD.N_Oxidizer = moles
        else:
            self.PD.N_Oxidizer = self.PD.phi.t / self.PD.phi.Value[i]
        self.PD.R_Oxidizer = SetSpecies(self)
    return self


def Define_I(self, species, i, **moles):
    self.PD.S_Inert = species
    if species:  # not empty, i.e., there is/are oxidizer/s in the mixture
        if moles:
            self.PD.N_Inert = moles
        else:
            self.PD.N_Inert = self.PD.phi.t / self.PD.phi.Value[i] * 79/21
        self.PD.R_Inert = SetSpecies(self)
    if self.PD.S_Oxidizer or self.PD.S_Inert:
        self.PS.strR_Oxidizer_and_Inert.append(ComputeProperties(
            self, self.PD.R_Oxidizer + self.PD.R_Inert, self.PD.pR.Value, self.PD.TR.Value))
    return self


def Define_FOI(self):
    R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert
    self.PS.strR.append(ComputeProperties(
        self, R, self.PD.pR.Value, self.PD.TR.Value))
    return self
