# -*- coding: utf-8 -*-
"""
CREATE STOICHIOMETRIC MATRIX

Created on Tue Jun 30 16:28:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np

def SetSpecies(self, S, N, T):
    M = self.C.M0.Value.copy()
    indexes = [self.S.NameSpecies.index(species) for species in S]
    R0 = self.C.R0

    for i, n in enumerate(N):
        hfi = self.strThProp[S[i]].hf / 1000.
        efi = self.strThProp[S[i]].ef / 1000.
        if len(self.strThProp[S[i]].T > 1):
            cPi = species_cP(S[i], T, self.strThProp)
            cVi = species_cV(S[i], T, self.strThProp)
            DeTi = species_DeT(S[i], T, self.strThProp)
            DhTi = species_DhT(S[i], T, self.strThProp)
            s0i = species_s0(S[i], T, self.strThProp)
            swtCondensed = self.strThProp[S[i]].swtCondensed
            mi = n * self.strThProp[S[i]].mm
            mmi = self.strThProp[S[i]].mm
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
            swtCondensed = self.strThProp[S[i]].swtCondensed
            mi = 0.
            mmi = self.strThProp[S[i]].mm
            pVi = 0.
        
        M[indexes[i], :] = np.concatenate(
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
