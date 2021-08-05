# Combustion PyToolbox
A Python based thermochemical code

🔜 `There is also a MATLAB version (coming soon)`

## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB and Python. The Python version solves, for the moment, six chemical equilibrium problems (`TP, HP, SP, TV, EV and SV transformations`; where T denotes temperature, P pressure, H enthalpy, S entropy, E internal energy and V volume), always assuming ideal gases in all cases.

---
⚠️ **NOTE**

At the moment, the Python version does not have all the capabilities that the MATLAB version has. I will continue with the development of this version adding all the capabilities that the MATLAB version has. I will also add a GUI using Qt and Pyside2.

---

The code computes the equilibrium composition by minimization of the Gibbs–Helmholtz free energy by using Lagrange multipliers, and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Results computed with **Combustion PyToolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program and CANTERA.

The `MATLAB version` `also solves` `incident and reflected planar shock waves`, as well as `ideal detonations` according to Chapman-Jouguet theory and overdriven detonations, assuming always ideal gases in all cases. Along with the plain code, the new tool has been `equipped with a Graphical User Interface` developed in MATLAB 2021 under AppDesigner.

This project is also part of my PhD.

## Installing

## Example
  ![](https://github.com/AlbertoCuadra/ThermochemicalCode_Python/blob/master/Validations/HP_CH4.svg)
  **Figure 1:** *An example of calculations given by Combustion PyToolbox: Adiabatic temperature and composition at constant pressure for methane with air.*
  
## Contributing

Please read [CONTRIBUTING.md](https://github.com/AlbertoCuadra/ThermochemicalCode_Python/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

We welcome your feedback, but do not guarantee support. Please send feedback or inquiries to [acuadra@ing.uc3m.es](mailto:acuadra@ing.uc3m.es)

Thank you for testing Combustion PyToolbox!

## Developers

* **[Alberto Cuadra-Lara](https://albertocuadra.netlify.app/)** - *Main Developer*
* **Marcos Vera** - *Developer*  

Grupo de Mecánica de Fluidos, Universidad Carlos III de Madrid, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.
