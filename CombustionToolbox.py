"""
COMBUSTION PYTOOLBOX

Type of problems:
    * TP ------> Equilibrium composition at defined T and p
    * HP ------> Adiabatic T and composition at constant p
    * SP ------> Isentropic compression/expansion to a specified p
    * TV ------> Equilibrium composition at defined T and constant v
    * EV ------> Adiabatic T and composition at constant v
    * SV ------> Isentropic compression/expansion to a specified v
    * SHOCK_I -> Planar incident shock wave
    * SHOCK_R -> Planar reflected shock wave
    * DET -----> Chapman-Jouget Detonation (CJ upper state)
    * DET_OVERDRIVEN -----> Overdriven Detonation    
    

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Thu Oct 1 13:10:00 2020
----------------------------------------------------------------------
"""
import time
import numpy as np
from Settings.Initialize import Initialize
from Settings.ListSpecies import ListSpecies
from Settings.Define_FOI import Define_F, Define_O, Define_I, Define_FOI
from Solver.Functions.Display.displayResults import displayResults
from Solver.Functions.Display.plotResults import plotResults
from Solver.Chemical_Equilibrium.SolveProblem import SolveProblem
from Solver.Functions.Transformation import set_transformation

def main():
    # LOAD DATABASES AND GLOBAL PARAMETERS
    app = Initialize()
    # LIST OF SPECIES
    app = ListSpecies(app, 'Soot formation')
    # app = ListSpecies(app, 'HC/02/N2 extended')
    # app = ListSpecies(app, 'Hydrogen')
    # app = ListSpecies(app, 'ideal_air')
    # PROBLEM TYPE AND CONDITIONS
    app.PD.ProblemType = 'HP' 
    
    set_transformation(app, 'TR', 300)  # [K]
    set_transformation(app, 'pR', 1.)   # [bar]
    set_transformation(app, 'TP', 2000) # [K]
    set_transformation(app, 'pP', 1.)   # [bar]
    
    app.PD.phi.Value = np.arange(0.5, 5, 0.01)  # [-]
    # COMPUTATIONS
    app.C.l_phi = len(app.PD.phi.Value)
    start = time.time()
    for i in range(app.C.l_phi):
        # DEFINE FUEL
        app = Define_F(app, {'CH4':1})
        # DEFINE OXIDIZER
        app = Define_O(app, {'O2':app.PD.phi.t / app.PD.phi.Value[i]})
        # DEFINE DILUENT/INERT
        app.PD.proportion_N2_O2 = 79/21
        app = Define_I(app, {'N2':app.PD.phi.t / app.PD.phi.Value[i] * app.PD.proportion_N2_O2})
        # COMPUTE PROPERTIES OF THE INITIAL MIXTURE
        app = Define_FOI(app, i)
        # SOLVE PROBLEM
        app = SolveProblem(app, i)
        # DISPLAY RESULTS
        displayResults(app, app.PS.strR[i], app.PS.strP[i])
    end = time.time()
    print('Execution time:', end - start, 'seconds')
    # PLOT RESULTS
    # display_species = ['CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb']
    display_species = app.S.LS
    plotResults(app, display_species=display_species, mintol=app.C.mintol_display)
    return app
    

if __name__ == '__main__':
    import cProfile, pstats
    profiler = cProfile.Profile()
    profiler.enable()
    # print(__doc__)
    main()
    profiler.disable()
    # stats = pstats.Stats(profiler).sort_stats('tottime')
    # stats.print_stats()
    # stats.dump_stats('/stats_file.dat')
    