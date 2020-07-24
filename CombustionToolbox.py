"""
COMBUSTION TOOLBOX

Type of problems:
    * TP ------> Equilibrium composition at defined T and p
    * HP ------> Adiabatic T and composition at constant p
    * SP ------> Isentropic compression/expansion to a specified p
    * TV ------> Equilibrium composition at defined T and constanxt v
    * EV ------> Adiabatic T and composition at constant v
    * SV ------> Isentropic compression/expansion to a specified v
    * SHOCK_I -> Planar incident shock wave
    * SHOCK_R -> Planar reflectet shock wave
    * DET -----> Chapman-Jouget Detonation (CJ upper state)
    * DET_OVERDRIVEN -----> Overdriven Detonation    
    

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
         
Last update Wen Jun 24 20:04:00 2020
----------------------------------------------------------------------
"""
import os
import time
from Settings.Initialize import Initialize
from Settings.Initialize_2 import Initialize_2
from Settings.MinorsProducts import MinorsProducts
from Settings.Define_FOI import *
from Solver.Functions.Display.displayResults import displayResults
from Solver.Chemical_Equilibrium.SolveProblemTP_TV import SolveProblemTP_TV

def main():
    # LOAD DATABASES AND GLOBAL PARAMETERS
    app = Initialize()
    # REACTION: COMPLETE OR INCOMPLETE
    #   Specify the type of the reaction: complete or incomplete (dissociation)
    app.PD.CompleteOrIncomplete = 'incomplete'  # incomplete (default)
    app.TN.factor_c = 1.0  # 1.0 (default)
    # MINORS PRODUCTS
    #   Specify the minority products to be considered in the product mixture (P) in
    #   addition to the major species (CO2, CO, H2O, H2, O2, N2, C(gr)).
    #   Moreover, He and Ar are always included in the database.
    #
    #   Specify in Setting/MinorsProducts.py
    #   Predefined:
    #       No minors products: ''
    #       * Hydrocarbons:               'HC/O2/N2 EXTENDED'
    #       * Soot formation:             'SOOT FORMATION'
    #       * Soot formation without CH4: 'SOOT FORMATION W/O CH4'
    #       * A bunch of minors products
    #         in case your are not sure
    #         which are possible minors
    #         products:                   'NASA ALL'
    #       * Air:                        'AIR'
    #       * Hydrogen:                   'HYDROGEN'
    #
    #   User definition:
    #       e.g., 'CH4, CO, O'

    app = MinorsProducts(app, 'Soot formation')
    # PROBLEM CONDITIONS

    # INITIALIZATION
    app = Initialize_2(app)
    # PROBLEM TYPE AND CONDITIONS
    app.PD.TR.Value = 300.  # [K]
    app.PD.pR.Value = 1.   # [bar]
    app.PD.phi.Value = [0.7]  # [-]
    
    app.PD.TP.Value = 2000
    # COMPUTATIONS
    app.C.l_phi = len(app.PD.phi.Value)
    for i in range(app.C.l_phi):
        # DEFINE FUEL
        app = Define_F(app, {'CH4':1})
        # DEFINE OXIDIZER
        app = Define_O(app, {'O2':app.PD.phi.t / app.PD.phi.Value[i]})
        # DEFINE DILUENT/INERT
        app = Define_I(app, {'N2':app.PD.phi.t / app.PD.phi.Value[i] * 79/21})
        # COMPUTE PROPERTIES OF THE INITIAL MIXTURE
        app = Define_FOI(app, i)
        # SOLVE PROBLEM
        # start = time.time()
        app.PS.strP.append(SolveProblemTP_TV(app, app.PS.strR[i], app.PD.phi.Value[i], app.PD.pR.Value, app.PD.TP.Value))
        # end = time.time()
        # DISPLAY RESULTS
        displayResults(app, app.PS.strR[i], app.PS.strP[i])
        # print(end - start)
    return app


if __name__ == '__main__':
    import cProfile, pstats
    profiler = cProfile.Profile()
    profiler.enable()
    # print(__doc__)
    # app = main()
    main()
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats('tottime')
    stats.print_stats()
    # stats.dump_stats('/stats_file.dat')
    
    