# Combustion PyToolbox
A Python based thermochemical code

## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB and Python. The Python version solves, for the moment, two chemical equilibrium problems (TP and HP transformations; where T denotes temperature, P pressure and H enthalpy), always assuming ideal gases in all cases. `The Python version does not have all the capabilities that the MATLAB version has, because it is much slower than this. When I or we (the community) fix this bottleneck, I will continue with the development of this version adding all the capabilities that the MATLAB version has. I will also add a GUI using Qt and Pyside2.`

The code computes the equilibrium composition by minimization of the Gibbs–Helmholtz free energy via Lagrange multipliers, and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Results computed with **Combustion PyToolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program and CANTERA.

This project is also part of my PhD. `The MATLAB version will be released soon.`  

The `MATLAB version` `solves` six chemical equilibrium problems (`TP, HP, SP, TV, EV and SV transformations`; where T denotes temperature, P pressure, H enthalpy, S entropy, E internal energy and V volume), `incident and reflected planar shock waves`, as well as `ideal detonations` according to Chapman-Jouguet theory and overdriven detonations, assuming always ideal gases in all cases. Along with the plain code, the new tool has been `equipped with a Graphical User Interface` (hereafter Combustion Toolbox) developed in MATLAB 2018-21 under AppDesigner.

## Installing

## Example
  ![Example computations given by Combustion PyToolbox: Adiabatic temperature and composition at constant pressure for methane.](https://github.com/AlbertoCuadra/ThermochemicalCode_Python/blob/master/Validations/HP_CH4.svg)
  
## Contributing

Please read [CONTRIBUTING.md](https://github.com/AlbertoCuadra/ThermochemicalCode_Python/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Developers

* **[Alberto Cuadra-Lara](https://albertocuadra.netlify.app/)** - *Main Developer*
* **Marcos Vera** - *Developer*  

Grupo de Mecánica de Fluidos, Universidad Carlos III, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.
