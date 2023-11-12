# -*- coding: utf-8 -*-
"""
Calculates the thermodynamic properties of any species included in the NASA database

Created on Wen Jun 24 20:04:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np
from utils.databases.set_element_matrix import set_element_matrix
from utils.databases.set_reference_form_of_elements_with_T_intervals import set_reference_form_of_elements_with_T_intervals
from utils.databases.isRefElm import isRefElm
from utils.databases.FullName2name import FullName2name
from utils.databases.detect_location_of_phase_specifier import detect_location_of_phase_specifier


def SpeciesThermProp(self, species, T, MassorMolar, echo):

    species = FullName2name(species)

    if species not in self.DB_master:
        if echo:
            print('Species %s does not exist as a field in DB_master structure' % species)
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

    name = self.DB_master[species].name
    FullName = self.DB_master[species].FullName
    comments = self.DB_master[species].comments
    ctTInt = self.DB_master[species].ctTInt
    txRefCode = self.DB_master[species].txRefCode
    txFormula = self.DB_master[species].txFormula
    swtCondensed = self.DB_master[species].swtCondensed
    mm = self.DB_master[species].mm
    Hf0 = self.DB_master[species].Hf0
    tRange = self.DB_master[species].tRange
    tExponents = self.DB_master[species].tExponents
    Hf298De10 = self.DB_master[species].Hf298De10

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
        a = self.DB_master[species].a
        b = self.DB_master[species].b
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

            DfG0 = H0 - T * S0
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
