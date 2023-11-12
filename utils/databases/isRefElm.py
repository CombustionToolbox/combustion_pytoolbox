# -*- coding: utf-8 -*-
"""
Checks if the species is a reference element (e.g., 'C(gr)')

Created on Wen Jun 24 20:04:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
from itertools import count


def isRefElm(reference, species, T):
    # Change lowercase 'l' to uppercase 'L' for Al, Cl, Tl, and Fl
    listL = ['Al', 'Cl', 'Tl', 'Fl']
    for specieL in listL:
        ind = species.find(specieL)
        if ind >= 0:
            species = species[:ind + 1] + species[ind +
                                                  1].replace('l', 'L') + species[ind + 2:]

    iRe = False
    REname = []
    """
    Look for entries in the Reference_form_of_elements_with_T_intervals list
    that partially match with the desired species and then check each one
    sucessivelly
    """
    ind = [i for i, j in zip(count(), reference) if species in j]
    for i in range(0, len(ind)):
        TentativeRefElm = reference[ind[i]]
        # Obtain temperature interval
        n1 = TentativeRefElm.find('[')
        n2 = TentativeRefElm.find('-')
        n3 = TentativeRefElm.find(']')
        T1 = float(TentativeRefElm[n1+1:n2])
        T2 = float(TentativeRefElm[n2+1:n3])
        if (T1 <= T) and (T <= T2):
            # Detect location of open parenthesis
            n_open_parenthesis = TentativeRefElm[0:n1-1].find('(')
            # Detect location of '2'
            n_two = TentativeRefElm[0:n1-1].find('2')
            # In case there is not a '2', the species is a eesentially a noble gas
            if n_open_parenthesis == -1 and n_two == -1:
                if TentativeRefElm[0:n1-1] == species:
                    iRe = True
                    REname = TentativeRefElm[0:n1-1]
                    break
            # If there are '2's the species may be in the reference state or
            # not (e.g. O2 is, but O is not)
            if n_two:
                if TentativeRefElm[0:n_two+1] == species:
                    iRe = True
                    REname = TentativeRefElm[0:n_two+1]
                    break
                if TentativeRefElm[0:n_two] == species:
                    REname = TentativeRefElm[0:n_two+1]
                    break
            # If there are opening parenthesis, the species is in condensed phase
            if n_open_parenthesis:
                if TentativeRefElm[0:n_open_parenthesis] == species:
                    iRe = True
                    REname = TentativeRefElm[0:n1-1]
                    break

    return [iRe, REname]
