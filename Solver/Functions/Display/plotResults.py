import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import (PchipInterpolator as pchip)

def plotMoles(self, x_vec, display_species, mintol, xlabel, ax=None, plt_args=[], plt_kwargs={}, **kwargs):
    if not ax:
        ax = plt.subplot()
    # Data
    struct = self.PS.strP
    Nx = len(x_vec)
    NS = len(display_species)
    y_matrix = np.zeros((NS, Nx))
    Ndisplay = set()
    # Display tolerance requirements
    for i in range(Nx):
        y_matrix[:, i] = struct[i].Xi
        for species in display_species:
            ind = self.S.LS.index(species)
            if struct[i].Xi[ind] > mintol:
                Ndisplay.add(ind)
    # Plot configuration
    plt.set_cmap('Spectral')
    ax.set_xlim(x_vec.min(), x_vec.max())
    ax.set_ylim(mintol, 1.0)
    ax.set_yscale('log')
    labels = np.array(display_species)
    labels = [labels[ndisplay] for ndisplay in Ndisplay]
    # Plot
    for j in Ndisplay:
        ax.plot(x_vec,y_matrix[j, :], *plt_args, **plt_kwargs, **kwargs)
    # Legend    
    ax.legend(labels=labels, loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.subplots_adjust(right=0.8)
    plt.show()
    return self

def plotResults(self, display_species=None, mintol=1e-16):
    phi = self.PD.phi.Value
    ProblemType = self.PD.ProblemType
    # Figure configuration
    plot_params = {'linewidth': 2}
    plot_args = ['-']
    if len(phi) > 1 and all(phi[1::] != phi[0]) and ProblemType != 'DET_OVERDRIVEN':
        plotMoles(self, phi, display_species, mintol, xlabel='Equivalence ratio', plt_args = plot_args, plt_kwargs = plot_params)
    