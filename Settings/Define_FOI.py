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


def Define_F(self, Species=None):
    if Species:  # not empty, i.e., there is/are fuel/s in the mixture
        self.PD.S_Fuel, self.PD.N_Fuel = unpack_dict(Species)
        self.PD.R_Fuel = SetSpecies(
            self, self.PD.S_Fuel, self.PD.N_Fuel, self.PD.TR.Value)
        self.PS.strR_Fuel = ComputeProperties(
            self, self.PD.R_Fuel, self.PD.pR.Value, self.PD.TR.Value)
        self.PD.Fuel.x = self.PS.strR_Fuel.NatomE[self.E.ind_C]
        self.PD.Fuel.y = self.PS.strR_Fuel.NatomE[self.E.ind_H]
        self.PD.Fuel.z = self.PS.strR_Fuel.NatomE[self.E.ind_O]
        self.PD.Fuel.w = self.PS.strR_Fuel.NatomE[self.E.ind_N]
        self.PD.phi.t = self.PD.Fuel.x + self.PD.Fuel.y / 4 - self.PD.Fuel.z / 2
    return self


def Define_O(self, Species=None):
    self.PD.S_Oxidizer, self.PD.N_Oxidizer = unpack_dict(Species)
    if Species:  # not empty, i.e., there is/are oxidizer/s in the mixture
        self.PD.R_Oxidizer = SetSpecies(
            self, self.PD.S_Oxidizer, self.PD.N_Oxidizer, self.PD.TR.Value)
    return self


def Define_I(self, Species=None):
    if Species:  # not empty, i.e., there is/are oxidizer/s in the mixture
        self.PD.S_Inert, self.PD.N_Inert = unpack_dict(Species)
        self.PD.R_Inert = SetSpecies(
            self, self.PD.S_Inert, self.PD.N_Inert, self.PD.TR.Value)
    else:
        self.PD.R_Inert = 0
    if self.PD.S_Oxidizer or self.PD.S_Inert:
        self.PS.strR_Oxidizer_and_Inert.append(ComputeProperties(
            self, self.PD.R_Oxidizer + self.PD.R_Inert, self.PD.pR.Value, self.PD.TR.Value))
    return self


def Define_FOI(self, i):
    R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert
    self.PS.strR.append(ComputeProperties(self, R, self.PD.pR.Value,
                                          self.PD.TR.Value))
    self.PS.strR[i].phi = self.PD.phi.Value[i]
    return self


def unpack_dict(dictionary):
    return (tuple(dictionary), list(dictionary.values()))
