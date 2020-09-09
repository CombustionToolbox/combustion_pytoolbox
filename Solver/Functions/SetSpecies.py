# -*- coding: utf-8 -*-
"""
CREATE STOICHIOMETRIC MATRIX

Created on Tue Jun 30 16:28:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np

def SetSpecies(self, Species, N, T):
    M = self.C.M0.Value.copy()
    R0 = self.C.R0

    for n, species in zip(N, Species):
        hfi = self.strThProp[species].hf / 1000.
        efi = self.strThProp[species].ef / 1000.
        if len(self.strThProp[species].T > 1):
            cPi = species_cP(species, T, self.strThProp)
            cVi = species_cV(species, T, self.strThProp)
            DeTi = species_DeT(species, T, self.strThProp)
            DhTi = species_DhT(species, T, self.strThProp)
            s0i = species_s0(species, T, self.strThProp)
            swtCondensed = self.strThProp[species].swtCondensed
            mi = n * self.strThProp[species].mm
            mmi = self.strThProp[species].mm
            if not swtCondensed:
                pVi = n * R0 * T / 100.  # For ideal gases
            else:
                pVi = 0.
        else:
            cPi = 0.
            cVi = 0.
            DeTi = 0.
            DhTi = 0.
            s0i = 0.
            swtCondensed = self.strThProp[species].swtCondensed
            mi = 0.
            mmi = self.strThProp[species].mm
            pVi = 0.
        
        M[self.S.List_Compute_Species.index(species), :] = np.concatenate(
            (n, n * np.array([hfi, DhTi, efi, DeTi, cPi, cVi, s0i]), pVi, swtCondensed, mi, mmi), axis=None)
    return M


def species_cP(species, T, strThProp):
    return strThProp[species].cPcurve(T)  # [J/mol K]


def species_cV(species, T, strThProp):
    return strThProp[species].cVcurve(T)  # [J/mol K]


def species_DeT(species, T, strThProp):
    return strThProp[species].DeTcurve(T) / 1000.  # [kJ/mol]


def species_DhT(species, T, strThProp):
    return strThProp[species].DhTcurve(T) / 1000.  # [kJ/mol]


def species_g0(species, T, strThProp):
    return strThProp[species].g0curve(T) / 1000.  # [kJ/mol]


def species_s0(species, T, strThProp):
    return strThProp[species].s0curve(T) / 1000.  # [kJ/mol K]

def equil_constant(DG0, TP, R0):
    return np.exp(-DG0 / (R0 * TP))