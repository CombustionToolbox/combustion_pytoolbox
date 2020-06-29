# -*- coding: utf-8 -*-
"""
LOAD/CALCULATE TABULATED DATA AND CONSTANTS

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import numpy as np
import re
from NASA_database.set_elements import *
from NASA_database.set_element_matrix import set_element_matrix
from NASA_database.ParseThermoInp import ParseThermoInp
from NASA_database.GenerateDatabase import GenerateDatabase

def Initialize():
    app = App()
    app.S.NameSpecies = tuple(app.strThProp.keys())
    app.S.NSpecies = len(app.S.NameSpecies)
    # Contained elements and species
    app.E = ContainedElements(app)
    app.E = Index_Evaluable_Elements(app.E)
    app.S = Index_Evaluable_Species(app.S)
    # Element Matrices
    app = Element_Matrices(app)
    # Guess initial calculation
    app.TN.guess = [2000., 2000., 0., 1.5, 2.]
    app.TN.ERRFT = 1e-5 # Tolerance SHOCK and Detonations numerical method
    app.TN.ERRFU = 1e-5 # Tolerance SHOCK and Detonations unmerical method
    
    return app

def ContainedElements(self):
    uniques = set() 
    for i, species in enumerate(self.S.NameSpecies):
        if any(species in element for element in self.E.Elements): # Case sensitive
            uniques.add(species)  
    
    self.E.Elements = tuple(uniques)
    self.E.ElementsUpper = UpperCase(self.E.Elements)
    self.E.NE = len(self.E.Elements)
    return self.E

def Index_Evaluable_Elements(self):
    self.ind_C = self.Elements.index('C')
    self.ind_H = self.Elements.index('H')
    self.ind_O = self.Elements.index('O')
    self.ind_N = self.Elements.index('N')
    self.ind_He = self.Elements.index('He')
    self.ind_Ar = self.Elements.index('Ar')
    
    return self

def Index_Evaluable_Species(self):
    self.ind_CO2 = self.NameSpecies.index('CO2')
    self.ind_CO = self.NameSpecies.index('CO')
    self.ind_H2O = self.NameSpecies.index('H2O')
    self.ind_H2 = self.NameSpecies.index('H2')
    self.ind_O2 = self.NameSpecies.index('O2')
    self.ind_N2 = self.NameSpecies.index('N2')
    self.ind_He = self.NameSpecies.index('He')
    self.ind_Ar = self.NameSpecies.index('Ar')
    self.ind_Cgr = self.NameSpecies.index('Cbgrb')
    
    self.List_fixed_Species = ['CO2','CO','H2O','H2','O2','N2','He','Ar','Cbgrb']
    self.idx_fixed = [self.ind_CO2, self.ind_CO, self.ind_H2O, self.ind_H2,\
        self.ind_O2, self.ind_N2, self.ind_He, self.ind_Ar, self.ind_Cgr]
    
    return self   

def Element_Matrices(self):
    self.C.A0.Value = np.zeros((self.S.NSpecies, self.E.NE))
    self.C.M0.Value = np.zeros((self.S.NSpecies, 12))
    for i, species in enumerate(self.S.NameSpecies):
        txFormula = self.strThProp[species].txFormula
        self.strThProp[species].Element_matrix = set_element_matrix(txFormula, self.E.ElementsUpper)
        ind_Elements, atoms = (self.strThProp[species].Element_matrix[0,:], self.strThProp[species].Element_matrix[1,:])
        for ind_Element, atom in zip(ind_Elements, atoms):
            self.C.A0.Value[i, int(ind_Element)] = atom
        self.C.M0.Value[i, 9] = self.strThProp[species].swtCondensed
    
    return self

# Definition of class App, which stores all the data necessary to execute the thermochemical code
class App:
    def __init__(self):
        self.E = self.Elements()
        self.S = self.Species()
        self.M = self.MinorsProducts()
        self.C = self.Constants()
        self.Misc = self.Miscelaneous()
        self.PD = self.ProblemDescription()
        self.PS = self.ProblemSolution()
        self.TN = self.TunningProperties()
        self.strMaster = ParseThermoInp(True) # False: complete DataBase; True: reduced DB
        self.strThProp = GenerateDatabase(self) # struct with tabulated data of selected species
    class Elements:
        def __init__(self):
            self.Description = "Data of the chemical elements"
            self.Elements, self.ElementsUpper, self.NE = set_elements()
    class Species:
        def __init__(self):
            self.Description = "Data of the chemical species"
            self.NameSpecies = []
            self.NSpecies = 0
    class MinorsProducts:
        def __init__(self):
            self.Description = "Data of minors products"
            self.display_species = []
            self.minor_products = []
    class Constants:
        def __init__(self):
            self.Description = "Constants and tolerances"
            self.R0 = 8.3144598 # [J/(K mol)]. Universal gas constant
            self.A0 = self.A0()
            self.M0 = self.M0()
            self.MassorMolar = 'mass'
            self.firstrow = True
            self.mintol_display = 1.e-8
            self.mintol = 1.e-5
            self.tolN = 1.e-14 # Tolerance for the segregated numerical method
            self.FLAG_FIRST = True 
        class A0:
            def __init__(self):
                self.Description = 'Molar Matrix: number of atoms of each element contained in each species'
                self.Value = []
        class M0:
            def __init__(self):
                self.Description = 'Matrix with properties of each species'
                self.Value = []
    class Miscelaneous:
        def __init__(self):
            self.Description = "Miscelaneous properties"
            self.Config = self.Config()
        class Config:
            def __init__(self):
                self.linewidth = 1.8
                self.fontsize = 18
    class ProblemDescription:
        def __init__(self):
            self.Description = "Problem description"
            self.CompleteOrIncomplete = "incomplete"            
    class ProblemSolution:
        def __init__(self):
            self.Description = "Problem solution"
    class TunningProperties:
        def __init__(self):
                self.Description = "Tunning properties"
                self.factor_c = 1.0
                
