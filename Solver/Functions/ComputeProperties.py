# -*- coding: utf-8 -*-
"""
DEFINE FUEL/OXIDIZER/INERT

Created on Wen Jun 30 12:49:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np
np.seterr(divide='ignore', invalid='ignore')


def ComputeProperties(self, SpeciesMatrix, p, T):
    struct = Struct(self, SpeciesMatrix, p, T)
    return struct


class Struct:
    def __init__(self, app, SpeciesMatrix, p, T):
        R0 = app.C.R0
        A0 = app.C.A0.Value
        ind_C, ind_H, ind_O, ind_N = [app.E.ind_C, app.E.ind_H, app.E.ind_O, app.E.ind_N]

        self.phi = None
        self.error_moles = None
        
        self.NatomE = sum(SpeciesMatrix[:, 0].reshape(
            SpeciesMatrix[:, 0].size, 1) * A0)
        self.x = self.NatomE[ind_C]
        self.y = self.NatomE[ind_H]
        self.z = self.NatomE[ind_O]
        self.w = self.NatomE[ind_N]
        self.N = sum(SpeciesMatrix[:, 0])
        self.hf = sum(SpeciesMatrix[:, 1])  # [kJ]
        self.DhT = sum(SpeciesMatrix[:, 2]) # [kJ]
        self.h = self.hf + self.DhT         # [kJ]
        self.ef = sum(SpeciesMatrix[:, 3])  # [kJ]
        self.ehT = sum(SpeciesMatrix[:, 4]) # [kJ]
        self.e = self.ef + self.ehT         # [kJ]
        self.S0 = sum(SpeciesMatrix[:, 7])  # [kJ/K]
        self.pv = sum(SpeciesMatrix[:, 8])  # [J/K]
        self.p = p                          # [bar]
        self.v = self.pv / self.p           # []
        self.T = T                          # [K]
        self.swtCond = SpeciesMatrix[:, 9]  # [-]
        self.mi = sum(SpeciesMatrix[:, 10]) * 1e-3            # [kg]
        self.rho = self.mi / self.v * 1e3                     # [kg/m3]
        self.Yi = SpeciesMatrix[:, 10] / self.mi * 1e-3       # [-]
        self.cP = sum(SpeciesMatrix[:, 5])  # [J/K]
        self.cV = sum(SpeciesMatrix[:, 6])  # [J/K]
        self.W = 1/np.nansum(self.Yi / SpeciesMatrix[:, 11])  # []
        Ni = SpeciesMatrix[:, 0]                              # [mol]
        self.Xi = Ni / self.N                                 # [-]
        posXi = [i for i, xi in enumerate(self.Xi) if xi > 0.]
        self.DSi = Ni[posXi] * \
            np.log(self.Xi[posXi] * p) * (1 - self.swtCond[posXi])
        self.DS = -R0 * sum(self.DSi)                         # [J/K]
        self.S = (self.S0 + self.DS * 1e-3) / self.mi         # [kJ/ kg-K]
        self.e = self.h - sum(Ni[posXi]) * R0 * T * 1e-3      # [kJ]
        self.gamma = self.cP / self.cV                        # [-]
        self.sound = np.sqrt(self.gamma * p * 1e5 / self.rho) # [m/s]
