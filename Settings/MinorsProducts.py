# -*- coding: utf-8 -*-
"""
Specify the minority products to be considered in the product mixture (P) in
addition to the major species (CO2, CO, H2O, H2, O2, N2, C(gr)).
Moreover, He and Ar are always included in the database.

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""


def MinorsProducts(self, type):
    if not type or type.upper() == 'NONE':
        self.M.minors_products = []
    elif type.upper() == 'HC/02/N2 EXTENDED':
        self.M.minors_products = ['OH', 'H', 'O', 'HO2', 'NO', 'HCO', 'CH4', 'CH3', 'HO2',
                                  'NO2', 'NH3', 'NH2', 'N', 'HCN', 'CN', 'N2O', 'C2', 'CH']
    elif type.upper() == 'SOOT FORMATION':
        self.M.minors_products = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2',
                                  'CH', 'CH2', 'CH3', 'CH4', 'CO', 'CO2', 'HCO', 'CH2OH', 'CH3O', 'CH3OH',
                                  'C2H', 'C2H4', 'C2H5', 'C2H6', 'HCCO', 'N', 'NH', 'NH2', 'NH3', 'NO',
                                  'NO2', 'N2O', 'HNO', 'CN', 'HCN', 'NCO', 'N2', 'Ar', 'C3H8', 'C2',
                                  'C2H2_acetylene', 'C6H6', 'C8H18_isooctane', 'C2H5OH', 'He', 'HNC', 'HNCO', 'NH2OH']
    elif type.upper() == 'SOOT FORMATION W/O CH4':
        self.M.minors_products = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2',
                                  'CH', 'CH2', 'CH3', 'CO', 'CO2', 'HCO', 'CH2OH', 'CH3O', 'CH3OH', 'C2H',
                                  'C2H4', 'C2H5', 'C2H6', 'HCCO', 'N', 'NH', 'NH2', 'NH3', 'NO', 'NO2',
                                  'N2O', 'HNO', 'CN', 'HCN', 'NCO', 'N2', 'Ar', 'C3H8', 'C2', 'C2H2_acetylene',
                                  'C6H6', 'C8H18_isooctane', 'C2H5OH', 'He', 'HNC', 'HNCO', 'NH2OH']
    elif type.upper() == 'NASA ALL':
        self.M.minors_products = ['CH3', 'CH4', 'CN', 'CO2', 'C2H', 'CH2CO_ketene',
                                  'C2H3_vinyl', 'C2H4', 'CH3COOH', 'C2H6', 'CH3OCH3', 'CNC', 'C2O',
                                  'C3H3_2_propynl', 'C3H4_cyclominus', 'C3H6_cyclominus', 'C3H6O_propanal',
                                  'C3H8', 'CNCOCN', 'C4H2_butadiyne', 'C4H6_1butyne', 'C4H8_1_butene',
                                  'C4H8_isobutene', 'C4H9_n_butyl', 'C4H9_t_butyl', 'C4N2', 'C5H8_cyclominus',
                                  'C5H11_pentyl', 'C5H12_i_pentane', 'C6H5_phenyl', 'C6H5OH_phenol',
                                  'C6H12_cyclominus', 'C7H7_benzyl', 'C7H14_1_heptene', 'C7H16_2_methylh',
                                  'C8H16_1_octene', 'C8H18_isooctane', 'C10H21_n_decyl', 'H', 'HCCN', 'HNCO',
                                  'HNO3', 'HCHO_formaldehy', 'H2O2', 'NCO', 'NH3', 'NO2', 'N2O', 'NH2NO2', 'N2O4',
                                  'N3H', 'O2', 'C2H5OHbLb', 'C6H6bLb', 'H2ObLb', 'CH', 'CH2OH', 'CH3OH',
                                  'CNN', 'COOH', 'C2H2_acetylene', 'ObCHb2O', 'CH3CN', 'C2H4O_ethylen_o',
                                  'OHCH2COOH', 'CH3N2CH3', 'CH3O2CH3', 'OCCN', 'C3', 'C3H4_allene', 'C3H5_allyl',
                                  'C3H6O_propylox', 'C3H7_n_propyl', 'C3H8O_1propanol', 'C3O2',
                                  'C4H4_1_3_cyclominus', 'C4H6_2butyne', 'C4H8_cis2_buten', 'C4H8_cyclominus',
                                  'C4H9_i_butyl', 'C4H10_n_butane', 'C5', 'C5H10_1_pentene', 'C5H11_t_pentyl',
                                  'CH3CbCH3b2CH3', 'C6H5O_phenoxy', 'C6H10_cyclominus', 'C6H13_n_hexyl', 'C7H8',
                                  'C7H15_n_heptyl', 'C8H8_styrene', 'C8H17_n_octyl', 'C9H19_n_nonyl', 'C12H9_o_bipheny',
                                  'HCN', 'HCCO', 'HNO', 'HO2', 'HCOOH', 'NH', 'NH2OH', 'NO3', 'NCN', 'N2H4',
                                  'N2O5', 'O', 'O3', 'N2H4bLb', 'CH2',
                                  'CH3O', 'CH3OOH', 'CO', 'C2', 'C2H2_vinylidene', 'HObCOb2OH', 'CH3CO_acetyl',
                                  'CH3CHO_ethanal', 'C2H5', 'C2H5OH', 'CCN', 'C2N2', 'C3H3_1_propynl', 'C3H4_propyne',
                                  'C3H6_propylene', 'C3H6O_acetone', 'C3H7_i_propyl', 'C3H8O_2propanol',
                                  'C4', 'C4H6_butadiene', 'C4H6_cyclominus', 'C4H8_tr2_butene',
                                  'C4H9_s_butyl', 'C4H10_isobutane', 'C5H6_1_3cyclominus', 'C5H10_cyclominus',
                                  'C5H12_n_pentane', 'C6H2', 'C6H6', 'C6H12_1_hexene', 'C6H14_n_hexane',
                                  'C7H8O_cresol_mx', 'C7H16_n_heptane', 'C8H10_ethylbenz', 'C8H18_n_octane',
                                  'C10H8_naphthale', 'C12H10_biphenyl', 'HCO', 'HNC', 'HNO2', 'H2', 'H2O',
                                  'N', 'NH2', 'NO', 'N2', 'N2H2', 'N2O3', 'N3', 'OH', 'CH3OHbLb', 'C6H5NH2bLb', 'He', 'Ar', 'C']
    elif type.upper() == 'AIR':
        self.M.minors_products = ['O2', 'N2', 'O', 'O3', 'N', 'NO', 'NO2', 'NO3', 'N2O3', 'N2O4', 'N3', 'C', 'CO', 'CO2',
                                  'Ar', 'CH4', 'CH3', 'CH', 'H2O', 'H2', 'H', 'He']
    elif type.upper() == 'HYDROGEN':
        self.M.minors_products = ['H', 'HNO', 'HNO3', 'H2O', 'NH', 'NH2OH', 'NO3', 'N2H2', 'N2O3', 'N3', 'OH', 'HNO2',
                                  'H2', 'N', 'NH3', 'NO2', 'N2O', 'N2H4', 'N2O5', 'O', 'O3', 'He', 'Ar', 'CO2', 'CO', 'O2', 'N2', 'HO2', 'NH2', 'H2O2',
                                  'N3H', 'NH2NO2']
    else:
        self.M.minors_products = type

    return self
