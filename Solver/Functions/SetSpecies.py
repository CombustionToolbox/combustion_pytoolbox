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
        hfi = self.strThProp[species].hf / 1000.
        efi = self.strThProp[species].ef / 1000.
        if len(self.strThProp[species].T > 1):
            tInterval = get_tInterval(species, T, self.strThProp)
            cPi = species_cP(species, T, self.strThProp, tInterval, R0)
            cVi = species_cV(species, T, self.strThProp, tInterval, R0)
            DeTi = species_DeT(species, T, self.strThProp, tInterval, R0)
            DhTi = species_DhT(species, T, self.strThProp, tInterval, R0)
            s0i = species_s0(species, T, self.strThProp, tInterval, R0)
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
        
        M[self.S.LS.index(species), :] = np.concatenate(
            (n, n * np.array([hfi, DhTi, efi, DeTi, cPi, cVi, s0i]), pVi, swtCondensed, mi, mmi), axis=None)
    return M

def get_tInterval(species, T, strThProp):
    # Select the appropriate temperature interval
    for i in range(0, strThProp[species].ctTInt):
        if (T >= strThProp[species].tRange[i][0]) and (T <= strThProp[species].tRange[i][1]):
            tInterval = i
    return tInterval

def species_cP(species, T, strThProp, tInterval, R0):
    return R0 * sum(strThProp[species].a[tInterval] * T**np.array(strThProp[species].tExponents[tInterval]))  # [J/mol K]


def species_cV(species, T, strThProp, tInterval, R0):
    cP = species_cP(species, T, strThProp, tInterval, R0)
    return cP - R0 # [J/mol K]


def species_DeT(species, T, strThProp, tInterval, R0):
    Tref = 298.15  # [K]
    H0 = species_DhT(species, T, strThProp, tInterval, R0) * 1000
    E0 = strThProp[species].ef + (H0 - strThProp[species].hf) - (1 - strThProp[species].swtCondensed) * R0 * (T - Tref)
    return (E0 - strThProp[species].ef) / 1000.  # [kJ/mol]


def species_DhT(species, T, strThProp, tInterval, R0):
    aux = np.array([-1, np.log(T), 1, 1/2, 1/3, 1/4, 1/5, 0])
    H0 = R0 * T * \
            (sum(strThProp[species].a[tInterval] * T**np.array(strThProp[species].tExponents[tInterval])
                 * aux) + strThProp[species].b[tInterval][0] / T)
    return (H0 - strThProp[species].hf) / 1000.  # [kJ/mol]


def species_g0(species, T, strThProp, tInterval, R0):
    DhT = species_DhT(species, T, strThProp, tInterval, R0) # [kJ/mol]
    H0 = DhT + strThProp[species].hf / 1000 # [kJ/mol]
    S0 = species_s0(species, T, strThProp, tInterval, R0) # [kJ/mol]
    return H0 - T * S0  # [kJ/mol]


def species_s0(species, T, strThProp, tInterval, R0):
    aux = np.array([-1/2, -1, np.log(T), 1, 1/2, 1/3, 1/4, 0])
    S0 = R0 * \
            (sum(strThProp[species].a[tInterval] * T**np.array(strThProp[species].tExponents[tInterval])
                 * aux) + strThProp[species].b[tInterval][1])
    return S0 / 1000.  # [kJ/mol K]


def equil_constant(DG0, TP, R0):
    return np.exp(-DG0 / (R0 * TP))