# -*- coding: utf-8 -*-
"""
Calculates the thermodynamic properties of any species included in the NASA database

Created on Wen Jun 24 20:04:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np
from NASA_database.set_element_matrix import set_element_matrix
from NASA_database.set_reference_form_of_elements_with_T_intervals import set_reference_form_of_elements_with_T_intervals
from NASA_database.isRefElm import isRefElm
from NASA_database.FullName2name import FullName2name
from NASA_database.detect_location_of_phase_specifier import detect_location_of_phase_specifier


def SpeciesThermProp(self, species, T, MassorMolar, echo):

    species = FullName2name(species)

    if species not in self.strMaster:
        if echo:
            print('Species %s does not exist as a field in strMaster structure' % species)
        txFormula = []
        mm = []
        Cp0 = []
        Cv0 = []
        Hf0 = []
        H0 = []
        Ef0 = []
        E0 = []
        S0 = []
        DfG0 = []

        return [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0]

    name = self.strMaster[species].name
    FullName = self.strMaster[species].FullName
    comments = self.strMaster[species].comments
    ctTInt = self.strMaster[species].ctTInt
    txRefCode = self.strMaster[species].txRefCode
    txFormula = self.strMaster[species].txFormula
    swtCondensed = self.strMaster[species].swtCondensed
    mm = self.strMaster[species].mm
    Hf0 = self.strMaster[species].Hf0
    tRange = self.strMaster[species].tRange
    tExponents = self.strMaster[species].tExponents
    Hf298De10 = self.strMaster[species].Hf298De10

    n_open_parenthesis = detect_location_of_phase_specifier(
        FullName)  # Detect the position of the phase specifier

    # Set Elements and reference form of elements with T intervals lists
    Element_matrix = set_element_matrix(txFormula, self.E.ElementsUpper)
    Reference_form_of_elements_with_T_intervals = set_reference_form_of_elements_with_T_intervals()

    """
    In order to compute the internal energy of formation from the enthalpy of
    formation of a given species, we must determine the change in moles of
    gases during the formation reaction of a mole of that species starting
    from the elements in their reference state. The only elements that are
    stable as diatomic gases are elements 1 (H), 7 (N), 8 (O), 9 (F), and 17
    (Cl). The remaining elements that are stable as (monoatomic) gases are
    the noble gases He (2), Ne (10), Ar (18), Kr (36), Xe (54), and Rn (86),
    which do not form any compound.
    """
    aux1 = np.array([[0], [6], [7], [8], [16]])  # Counting 0
    aux2 = np.array([[1], [9], [17], [35], [53], [87]])  # Counting 0
    Delta_n_per_mole = sum(Element_matrix[0, :] == aux1)/2 +\
        sum(Element_matrix[0, :] == aux2)
    Delta_n = 1. - swtCondensed - \
        np.dot(Delta_n_per_mole, Element_matrix[1, :])

    R0 = self.C.R0
    """
    Check if there is at least one temperature interval and, in that case,
    check that the specified temperature is within limits. If it is not, then
    abort, otherwise keep on running
    """
    if ctTInt > 0:
        a = self.strMaster[species].a
        b = self.strMaster[species].b
        Tref = 298.15  # [K]

        if (T < tRange[0][0]) or (T > tRange[ctTInt-1][1]) and echo:
            print('T out of range [%.2f - %.2f] [K] for %s' %
                  (tRange[0][0], tRange[ctTInt][1], FullName))
            Cp0 = []
            Cv0 = []
            H0 = []
            Ef0 = []
            E0 = []
            S0 = []
            DfG0 = []
            return [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0]

        # Select the appropriate temperature interval
        for i in range(0, ctTInt):
            if (T >= tRange[i][0]) and (T <= tRange[i][1]):
                tInterval = i

        """ 
        Compute the thermochemical data at the specified temperature using
        the polynomial coefficients in the selected temperature interval. All
        magnitudes are computed in a per mole basis  
        """
        Cp0 = R0 * sum(a[tInterval] * T**np.array(tExponents[tInterval]))
        Cv0 = Cp0 - R0
        aux = np.array([-1, np.log(T), 1, 1/2, 1/3, 1/4, 1/5, 0])
        H0 = R0 * T * \
            (sum(a[tInterval] * T**np.array(tExponents[tInterval])
                 * aux) + b[tInterval][0] / T)
        Ef0 = Hf0 - Delta_n * R0 * Tref
        E0 = Ef0 + (H0 - Hf0) - (1 - swtCondensed) * R0 * (T - Tref)
        aux = np.array([-1/2, -1, np.log(T), 1, 1/2, 1/3, 1/4, 0])
        S0 = R0 * \
            (sum(a[tInterval] * T**np.array(tExponents[tInterval])
                 * aux) + b[tInterval][1])

        """
        Compute the standar gibbs free energy of formation at the specified
        temperature. This enforces us to consider explicitely the formation
        reaction from the elements in their reference states at room
        temperature, unless the species is precisely an element in its
        reference state, in which case the standard gibbs free energy of
        formation is identically zero.
        """
        [iRe, REname] = isRefElm(
            Reference_form_of_elements_with_T_intervals, FullName[0:n_open_parenthesis], T)
        if not iRe:
            if echo:
                print(f'{FullName} is not Ref-Elm')

            GP = H0 - T * S0
            GR = np.zeros((1, len(Element_matrix[1, :])))
            for i in range(0, len(Element_matrix[1, :])):
                nu_i = Element_matrix[1, i]
                [iRe_i, REname_i] = isRefElm(
                    Reference_form_of_elements_with_T_intervals, self.E.ElementsUpper[int(Element_matrix[0, i])], T)
                [txFormula_i, mm_i, Cp0_i, Cv0_i, Hf0_i, H0_i, Ef0_i, E0_i,
                    S0_i, DfG0_i] = SpeciesThermProp(self, REname_i, T, 'molar', 0)
                GR[0, i] = nu_i * (H0_i - T * S0_i)
                if any(Element_matrix[0, i] == [0, 6, 7, 8, 16, 34]):
                    GR[0, i] = GR[0, i]/2.

            GR = GR.sum()
            DfG0 = GP - GR
        else:
            if echo:
                print(f'{REname} is Ref-Elm')
            DfG0 = 0.

        if MassorMolar == 'mass':
            Cp0 = Cp0 / (mm / 1000)
            Cv0 = Cv0 / (mm / 1000)
            Hf0 = Hf0 / (mm / 1000)
            Ef0 = Ef0 / (mm / 1000)
            H0 = H0 / (mm / 1000)
            S0 = S0 / (mm / 1000)
            if not swtCondensed:
                DfG0 = DfG0 / (mm / 1000)
            else:
                DfG0 = []

    else:
        """        
        If the species is only a reactant determine it's reference temperature
        Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
        K. For cryogenic liquids, assigned enthalpies are given at their boiling
        points instead of 298.15 K
        """
        if T != tRange[0]:
            print('T out of range [%.2f - %.2f] [K] for %s' %
                  (tRange[0][0], tRange[ctTInt][1], FullName))
            Cp0 = []
            Cv0 = []
            H0 = []
            Ef0 = Hf0 - Delta_n * R0 * tRange[0]
            E0 = []
            S0 = []
            DfG0 = []
        else:
            Tref = tRange[0]
            Cp0 = 0.
            Cv0 = 0.
            H0 = 0.
            E0 = 0.
            Ef0 = Hf0 - Delta_n * R0 * Tref

        if MassorMolar == 'mass':
            Hf0 = Hf0 / (mm / 1000)
            Ef0 = Ef0 / (mm / 1000)

    return [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0]
