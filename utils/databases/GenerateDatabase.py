# -*- coding: utf-8 -*-
"""
Generate tabulated database (NASA)

Created on Wen Jun 24 20:04:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np
import pickle
from utils.databases.FullName2name import FullName2name
from utils.databases.SpeciesThermProp import SpeciesThermProp
from time import time


def GenerateDatabase(self):  # self is a dictionary with Master Database
    filepath = 'Databases/DB.pkl'
    if not existFile(filepath):
        print('Generating tabulated NASA database ...', end=" ")

        Species = ['C', 'CH3', 'CH4', 'CN', 'CO2', 'C2H', 'CH2CO_ketene',
                   'C2H3_vinyl', 'C2H4', 'CH3COOH', 'C2H6', 'CH3OCH3', 'CNC', 'C2O',
                   'C3H3_2_propynl', 'C3H6O_propanal', 'C3H8', 'CNCOCN', 'C4H2_butadiyne',
                   'C4H6_1butyne', 'C4H8_1_butene', 'C4H8_isobutene', 'C4H9_n_butyl',
                   'C4H9_t_butyl', 'C4N2', 'C5H11_pentyl', 'C5H12_i_pentane', 'C6H5_phenyl',
                   'C6H5OH_phenol', 'C7H7_benzyl', 'C7H14_1_heptene', 'C7H16_2_methylh',
                   'C8H16_1_octene', 'C8H18_isooctane', 'C10H21_n_decyl', 'H', 'HCCN',
                   'HNCO', 'HNO3', 'HCHO_formaldehy', 'H2O2', 'NCO', 'NH3', 'NO2', 'N2O',
                   'NH2NO2', 'N2O4', 'N3H', 'O2', 'Cbgrb', 'C2H5OHbLb', 'C6H6bLb', 'H2ObLb',
                   'CH', 'CH2OH', 'CH3OH', 'CNN', 'COOH', 'C2H2_acetylene', 'ObCHb2O', 'CH3CN',
                   'C2H4O_ethylen_o', 'OHCH2COOH', 'CH3N2CH3', 'CH3O2CH3', 'OCCN', 'C3',
                   'C3H4_allene', 'C3H5_allyl', 'C3H6O_propylox', 'C3H7_n_propyl',
                   'C3H8O_1propanol', 'C3O2', 'C4H6_2butyne', 'C4H8_cis2_buten',
                   'C4H9_i_butyl', 'C4H10_n_butane', 'C5', 'C5H10_1_pentene', 'C5H11_t_pentyl',
                   'CH3CbCH3b2CH3', 'C6H5O_phenoxy', 'C6H13_n_hexyl', 'C7H8',
                   'C7H15_n_heptyl', 'C8H8_styrene', 'C8H17_n_octyl', 'C9H19_n_nonyl',
                   'C12H9_o_bipheny', 'HCN', 'HCCO', 'HNO', 'HO2', 'HCOOH', 'NH', 'NH2OH',
                   'NO3', 'NCN', 'N2H4', 'N2O5', 'O', 'O3', 'N2H4bLb', 'CH2', 'CH3O', 'CH3OOH',
                   'CO', 'C2', 'C2H2_vinylidene', 'HObCOb2OH', 'CH3CO_acetyl',
                   'CH3CHO_ethanal', 'C2H5', 'C2H5OH', 'CCN', 'C2N2', 'C3H3_1_propynl',
                   'C3H4_propyne', 'C3H6_propylene', 'C3H6O_acetone', 'C3H7_i_propyl',
                   'C3H8O_2propanol', 'C4', 'C4H6_butadiene', 'C4H8_tr2_butene',
                   'C4H9_s_butyl', 'C4H10_isobutane', 'C5H12_n_pentane', 'C6H2', 'C6H6',
                   'C6H12_1_hexene', 'C6H14_n_hexane', 'C7H8O_cresol_mx', 'C7H16_n_heptane',
                   'C8H10_ethylbenz', 'C8H18_n_octane', 'C10H8_naphthale', 'C12H10_biphenyl',
                   'HCO', 'HNC', 'HNO2', 'H2', 'H2O', 'N', 'NH2', 'NO', 'N2', 'N2H2', 'N2O3', 'N3',
                   'OH', 'CH3OHbLb', 'C6H5NH2bLb', 'He', 'Ar', 'Cbgrb', 'F2', 'F']

        DB = {}
        for ind, FullSpecies in enumerate(Species):
            aux = DataBase()
            species = FullName2name(FullSpecies)
            if species in self.DB_master:
                ctTInt = self.DB_master[species].ctTInt
                tRange = self.DB_master[species].tRange
                swtCondensed = self.DB_master[species].swtCondensed
                if ctTInt > 0:
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(
                        self, species, 298.15, 'molar', 0)
                    aux.name = species
                    aux.FullName = FullSpecies
                    aux.txFormula = txFormula
                    aux.mm = mm
                    aux.hf = Hf0
                    aux.ef = Ef0
                    aux.swtCondensed = swtCondensed

                    T_vector = np.array([])
                    DhT_vector = np.array([])
                    DeT_vector = np.array([])
                    s0_vector = np.array([])
                    cp_vector = np.array([])
                    cv_vector = np.array([])
                    g0_vector = np.array([])

                    Tmin = max(tRange[0][0], 200.)
                    Tmax = min(tRange[ctTInt-1][1], 6000.)

                    for T in np.concatenate((np.linspace(Tmin, 298.15, 4),  np.linspace(350., Tmax, 60))):
                        [_, _, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(
                            self, FullSpecies, T, 'molar', 0)
                        T_vector = np.hstack((T_vector, T))
                        DhT_vector = np.hstack([DhT_vector, H0 - Hf0])
                        DeT_vector = np.hstack([DeT_vector, E0 - Ef0])
                        s0_vector = np.hstack([s0_vector, S0])
                        cp_vector = np.hstack([cp_vector, Cp0])
                        cv_vector = np.hstack([cv_vector, Cv0])
                        g0_vector = np.hstack([g0_vector, DfG0])

                    aux.T = T_vector
                    aux.DhT = DhT_vector
                    aux.DeT = DeT_vector
                    aux.s0 = s0_vector
                    aux.cp = cp_vector
                    aux.cv = cv_vector
                    aux.g0 = g0_vector
                    aux.a = self.DB_master[species].a
                    aux.b = self.DB_master[species].b
                    aux.tExponents = self.DB_master[species].tExponents
                    aux.ctTInt = ctTInt
                    aux.tRange = tRange
                else:
                    Tref = aux.tRange[0]

                    aux.name = species
                    aux.FullName = FullSpecies
                    aux.txFormula = txFormula
                    aux.mm = mm
                    aux.hf = Hf0
                    aux.ef = Ef0
                    aux.swtCondensed = swtCondensed
                    aux.T = Tref

                DB.update({aux.name: aux})
            else:
                print(f'Species {FullSpecies} does not exist as a field in strMaster structure')
        # Save StrThprop
        f = open("Databases/DB.pkl", "wb")
        pickle.dump(DB, f)
        f.close()

        print('OK!')
    else:  # Load StrMaster reduced
        with open(filepath, 'rb') as f:
            DB = pickle.load(f)
        print(f'NASA tabulated database loaded from {filepath}')

    return DB


def existFile(filepath):
    import os
    return os.path.exists(filepath)


class DataBase():
    def __init__(self):
        self.FullName = []
        self.name = []
        self.txFormula = []
        self.swtCondensed = []
        self.mm = []
        self.hf = []
        self.ef = []
        self.T = []
        self.DhT = []
        self.DeT = []
        self.s0 = []
        self.s = []
        self.cp = []
        self.cv = []
        self.g0 = []
        self.a = []
        self.b = []
        self.tRange = []
        self.tExponents = []
