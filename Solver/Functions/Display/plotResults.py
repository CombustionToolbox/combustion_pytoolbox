import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Paired_12
from Solver.Functions.struct2vector import struct2vector

def plotMoles(self, x_vec, display_species, mintol, ax=None, plt_args=[], plt_kwargs={}, **kwargs):
    plt.ioff()
    if not ax:
        self.Misc.CounterPlots += 1 
        ax = plt.subplot()
    # Data
    struct = self.PS.strP
    Nx = len(x_vec)
    NS = len(struct[0].Xi)
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
    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9, hspace=0.5)
    title  = kwargs.pop('title')
    xlabel = kwargs.pop('xlabel')
    ylabel = kwargs.pop('ylabel')
    
    Ldisplay = len(Ndisplay)
    maxLdisplay = 12
    if Ldisplay > maxLdisplay:
        NUM_COLORS = maxLdisplay
    else:
        NUM_COLORS = Ldisplay

    LINE_STYLES = ['solid', 'dashed', 'dotted']
    NUM_STYLES = len(LINE_STYLES)

    cmap = matplotlib.colors.ListedColormap(Paired_12.mpl_colors)
    colors = cmap(np.linspace(0, 1.0, NUM_COLORS))
    # plt.set_cmap(cmap)
    # sns.reset_orig()  # get default matplotlib styles back
    # colors = sns.color_palette('husl', n_colors=NUM_COLORS)  # a list of RGB tuples
    
    ax.set_xlim(x_vec.min(), x_vec.max())
    ax.set_ylim(mintol, 1.0)
    ax.set_yscale('log')
    ax.set_title(title)
    ax.xaxis.set_label_text(xlabel)
    ax.yaxis.set_label_text(ylabel)
    labels = np.array(display_species)
    labels = [labels[ndisplay] for ndisplay in Ndisplay]
    # Plot
    k = 0
    z = 0
    for j in Ndisplay:
        lines = ax.plot(x_vec, y_matrix[j, :], *plt_args, **plt_kwargs, **kwargs)
        if Ldisplay > maxLdisplay:
            lines[0].set_color(colors[k])
            lines[0].set_linestyle(LINE_STYLES[z%NUM_STYLES])
            k = k + 1
            if k == 12:
                k = 0
                z = z + 1

    # Legend    
    ax.legend(labels=labels, loc='upper left', bbox_to_anchor=(1.05, 1))

    return self

def plotFigure(self, x_vec, y_vec, ax=None, plt_args=[], plt_kwargs={}, **kwargs):
    plt.ioff()
    if not ax:
        self.Misc.CounterPlots += 1 
        plt.figure(self.Misc.CounterPlots)
        ax = plt.gca()
    # Plot configuration
    title  = kwargs.pop('title')
    xlabel = kwargs.pop('xlabel')
    ylabel = kwargs.pop('ylabel')
    ax.set_xlim(x_vec.min(), x_vec.max())
    ax.set_ylim(y_vec.min(), 1.05*y_vec.max())
    ax.set_title(title)
    ax.xaxis.set_label_text(xlabel)
    ax.yaxis.set_label_text(ylabel)
    # Plot
    ax.plot(x_vec, y_vec, *plt_args, **plt_kwargs, **kwargs)
    
    return self

def plotResults(self, display_species=None, mintol=1e-16):
    phi = self.PD.phi.Value
    ProblemType = self.PD.ProblemType
    # Figure configuration
    plot_args = ['-']
    plot_params = {'linewidth': 2} # Parameters for plot lines
    if len(phi) > 1 and all(phi[1::] != phi[0]) and ProblemType != 'DET_OVERDRIVEN':
        plot_conf = {'title': 'Molar fraction', 'xlabel':r'$\phi$', 'ylabel': r'$X_i$'}        
        plotMoles(self, phi, display_species, mintol, plt_args = plot_args, plt_kwargs = plot_params, **plot_conf)
        if self.PD.ProblemType.upper() != ('TP' or 'TV'):
            plot_conf = {'title': 'Adiabatic temperature', 'xlabel':r'$\phi$', 'ylabel': r'$T [K]$'}
            plotFigure(self, phi, struct2vector(self.PS.strP, 'T'), plt_args = plot_args, plt_kwargs = plot_params, **plot_conf)
    
    plt.show()
