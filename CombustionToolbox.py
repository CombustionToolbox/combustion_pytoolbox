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
"""
import os
from Settings.Initialize import App

def main():
    app = App()
    return app
    
if __name__ == '__main__':
    app = main()
