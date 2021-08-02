# Combustion PyToolbox
A Python based thermochemical code

## Introduction
As a first step towards the development of a wider-scope thermochemical tool, in this work we present a thermochemical code with application to gaseous combustion problems recently implemented by the authors in MATLAB and Python. The Python version solves, for the moment, a chemical equilibrium problem (TP transformations; where T denotes temperature, P pressure), always assuming ideal gases. `The Python version does not have all the capabilities that the MATLAB version has, because it is much slower than this. When I or we (the community) fix this bottleneck, I will continue with the development of this version adding all the capabilities that the MATLAB version has. I will also add a GUI using Qt5 and Pyside2.`

The code can compute the equilibrium composition by minimization of the Gibbs–Helmholtz free energy or using a using equilibrium constants (segregated method), and employs NASA’s 9-coefficient polynomial fits to evaluate the thermodynamic properties. Results computed with **Combustion-Toolbox** have been validated against, and are in good agreement with, NASA’s Chemical Equilibrium with Applications (CEA) program and CANTERA.

## Installing

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **[Alberto Cuadra-Lara](https://albertocuadra.netlify.app/)** - *Main Developer*
* **Marcos Vera** - *Developer*  

Grupo de Mecánica de Fluidos, Universidad Carlos III, Av. Universidad 30, 28911, Leganés, Spain

See also the list of [contributors](https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md) who participated in this project.
