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
from NASA_database.set_elements import set_elements
from NASA_database.ParseThermoInp import ParseThermoInp
from NASA_database.GenerateDatabase import GenerateDatabase

def Initialize():
    app = App()
    app.S.NameSpecies = app.strThProp.keys()
    app.S.NSpecies = len(app.S.NameSpecies)
    # Contained elements
    return app

def ContainedElements(self):
    for i, species in enumerate(self.S.NameSpecies):
        species = re.sub('plus', '')
        species = re.sub('minus', '')
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
            self.Elements, self.NE = set_elements()
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
                
