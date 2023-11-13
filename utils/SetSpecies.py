# -*- coding: utf-8 -*-
"""
CREATE STOICHIOMETRIC MATRIX

Created on Tue Jun 30 16:28:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
"""
import numpy as np

def SetSpecies(self, Species, N, T):
    M = self.C.M0.Value.copy()
    R0 = self.C.R0

    for n, species in zip(N, Species):
        hfi = self.DB[species].hf / 1000.
        efi = self.DB[species].ef / 1000.
        if len(self.DB[species].T > 1):
            tInterval = get_tInterval(species, T, self.DB)
            cPi = species_cP(species, T, self.DB, tInterval, R0)
            cVi = species_cV(species, T, self.DB, tInterval, R0)
            DeTi = species_DeT(species, T, self.DB, tInterval, R0)
            DhTi = species_DhT(species, T, self.DB, tInterval, R0)
            s0i = species_s0(species, T, self.DB, tInterval, R0)
            swtCondensed = self.DB[species].swtCondensed
            mi = n * self.DB[species].mm
            mmi = self.DB[species].mm
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
            swtCondensed = self.DB[species].swtCondensed
            mi = 0.
            mmi = self.DB[species].mm
            pVi = 0.
        
        M[self.S.LS.index(species), :] = np.concatenate(
            (n, n * np.array([hfi, DhTi, efi, DeTi, cPi, cVi, s0i]), pVi, swtCondensed, mi, mmi), axis=None)
    return M

def get_tInterval(species, T, DB):
    # Select the appropriate temperature interval
    for i in range(0, DB[species].ctTInt):
        if (T >= DB[species].tRange[i][0]) and (T <= DB[species].tRange[i][1]):
            tInterval = i
    return tInterval

def species_cP(species, T, DB, tInterval, R0):
    return R0 * sum(DB[species].a[tInterval] * T**np.array(DB[species].tExponents[tInterval]))  # [J/mol K]


def species_cV(species, T, DB, tInterval, R0):
    cP = species_cP(species, T, DB, tInterval, R0)
    return cP - R0 # [J/mol K]


def species_DeT(species, T, DB, tInterval, R0):
    Tref = 298.15  # [K]
    H0 = species_DhT(species, T, DB, tInterval, R0) * 1000
    E0 = DB[species].ef + (H0 - DB[species].hf) - (1 - DB[species].swtCondensed) * R0 * (T - Tref)
    return (E0 - DB[species].ef) / 1000.  # [kJ/mol]


def species_DhT(species, T, DB, tInterval, R0):
    aux = np.array([-1, np.log(T), 1, 1/2, 1/3, 1/4, 1/5, 0])
    H0 = R0 * T * \
            (sum(DB[species].a[tInterval] * T**np.array(DB[species].tExponents[tInterval])
                 * aux) + DB[species].b[tInterval][0] / T)
    return (H0 - DB[species].hf) / 1000.  # [kJ/mol]


def species_g0(species, T, DB, tInterval, R0):
    DhT = species_DhT(species, T, DB, tInterval, R0) # [kJ/mol]
    H0 = DhT + DB[species].hf / 1000 # [kJ/mol]
    S0 = species_s0(species, T, DB, tInterval, R0) # [kJ/mol]
    return H0 - T * S0  # [kJ/mol]


def species_s0(species, T, DB, tInterval, R0):
    aux = np.array([-1/2, -1, np.log(T), 1, 1/2, 1/3, 1/4, 0])
    S0 = R0 * \
            (sum(DB[species].a[tInterval] * T**np.array(DB[species].tExponents[tInterval])
                 * aux) + DB[species].b[tInterval][1])
    return S0 / 1000.  # [kJ/mol K]


def equil_constant(DG0, TP, R0):
    return np.exp(-DG0 / (R0 * TP))