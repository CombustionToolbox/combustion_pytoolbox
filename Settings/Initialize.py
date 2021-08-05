# -*- coding: utf-8 -*-
"""
LOAD/CALCULATE TABULATED DATA AND CONSTANTS

Created on Mon Jun 22 12:15:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
"""

from NASA_database.set_elements import *
from NASA_database.ParseThermoInp import ParseThermoInp
from NASA_database.GenerateDatabase import GenerateDatabase


def Initialize():
    app = App()
    app.S.NameSpecies = tuple(app.strThProp.keys())
    app.S.NSpecies = len(app.S.NameSpecies)
    # Contained elements
    app.E = ContainedElements(app)
    app.E = Index_Evaluable_Elements(app.E)
    # Guess initial calculation
    app.TN.guess = [2000., 2000., 0., 1.5, 2.]
    app.TN.ERRFT = 1e-5  # Tolerance SHOCK and Detonations numerical method
    app.TN.ERRFU = 1e-5  # Tolerance SHOCK and Detonations numerical method

    return app


def ContainedElements(self):
    uniques = set()
    for species in self.S.NameSpecies:
        if any(species in element for element in self.E.Elements):  # Case sensitive
            uniques.add(species)

    self.E.Elements = sorted(uniques)
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

# Definition of class App, which stores all the data necessary to execute the thermochemical code


class App:
    def __init__(self):
        self.E = self.Elements()
        self.S = self.Species()
        self.C = self.Constants()
        self.Misc = self.Miscelaneous()
        self.PD = self.ProblemDescription()
        self.PS = self.ProblemSolution()
        self.TN = self.TunningProperties()
        # False: complete DataBase; True: reduced DB
        self.strMaster = ParseThermoInp(False)
        # struct with tabulated data of selected species
        self.strThProp = GenerateDatabase(self)

    class Elements:
        def __init__(self):
            self.Description = "Data of the chemical elements"
            self.Elements, self.ElementsUpper, self.NE = set_elements()

    class Species:
        def __init__(self):
            self.Description = "Data of the chemical species"
            self.NameSpecies = []
            self.NSpecies = 0
            self.NG = 0
            self.NS = 0

    class Constants:
        def __init__(self):
            self.Description = "Constants and tolerances"
            self.R0 = 8.3144598  # [J/(K mol)]. Universal gas constant
            self.A0 = self.A0()
            self.M0 = self.M0()
            self.MassorMolar = 'mass'
            self.firstrow = True
            self.mintol_display = 1e-10
            self.mintol = 1e-5
            self.tolN = 1e-10 # Tolerance of the segregated numerical method
            self.tolPhiSoot = 1e-6  # Tolerance of the soot formation equivalence ratio numerical method

        class A0:
            def __init__(self):
                self.Description = 'Stoichiometric Matrix: number of atoms of each element contained in each species'
                self.Value = None

        class M0:
            def __init__(self):
                self.Description = 'Matrix with properties of each species'
                self.Value = None
                
        class N0:
            def __init__(self):
                self.Description = 'Reduced Matrix with number of moles and swtCondensated of each species'
                self.Value = None

    class Miscelaneous:
        def __init__(self):
            self.Description = "Miscellaneous properties"
            self.Config = self.Config()
            self.FLAG_FIRST = True
            self.CounterPlots = 0
        class Config:
            def __init__(self):
                self.linewidth = 1.8
                self.fontsize = 18

    class ProblemDescription:
        def __init__(self):
            self.Description = "Problem description"
            self.CompleteOrIncomplete = "incomplete"
            self.ProblemType = None
            self.R_Fuel = None
            self.R_Oxidizer = None
            self.R_Inert = None
            self.phi = self.Phi()
            self.Fuel = self.Fuel()
            self.TR = self.TR()
            self.TP = self.TP()
            self.pR = self.pR()
            self.pP = self.pP()
            self.vP_vR = self.vP_vR()
            self.u1 = self.u1()
            self.overdriven = self.overdriven()
            self.proportion_N2_O2 = None
        class Phi:
            def __init__(self):
                self.Description = "Equivalence ratio"
                self.Value = 1.0
                self.t = 1.0  # "Theoretical value: phi = phi_t / phi.st"

        class Fuel:
            def __init__(self):
                self.x = None  # C atoms in the fuel mixture
                self.y = None  # H atoms in the fuel mixture
                self.z = None  # O atoms in the fuel mixture
                self.w = None  # N atoms in the fuel mixture
                self.eps = 0.

        class TR:
            def __init__(self):
                self.Description = 'Temperature of reactants'
                self.Value = None
        
        class TP:
            def __init__(self):
                self.Description = 'Temperature of products'
                self.Value = None

        class pR:
            def __init__(self):
                self.Description = 'Pressure of reactants'
                self.Value = None
        class pP:
            def __init__(self):
                self.Description = 'Pressure of products'
                self.Value = None
        class vP_vR:
            def __init__(self):
                self.Description = 'Volume relation Products/Reactants'
                self.Value = None
        class u1:
            def __init__(self):
                self.Description = 'Incident shock velocity'
                self.Value = None        
        class overdriven:
            def __init__(self):
                self.Description = 'Overdriven shock velocity'
                self.Value = None        
    
    class ProblemSolution:
        def __init__(self):
            self.Description = "Problem solution"
            self.strR_Fuel = []  # Computations of the reactant fuel
            self.strR_Oxidizer_and_Inert = []  # Computations of the reactant fuel
            self.strR = []  # Computations of the reactants
            self.strP = []  # Computations of the products

    class TunningProperties:
        def __init__(self):
            self.Description = "Tunning properties"
            self.factor_c = 1.0
            self.ERRFT = 1e-6
            self.ERRFV = 1e-6
            self.ERRFU = 1e-4
