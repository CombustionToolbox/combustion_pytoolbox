"""
COMPUTE CHEMICAL EQUILIBRIUM USING THE GENERALIZED GIBBS MINIMIZATION METHOD
   
@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D22, Universidad Carlos III de Madrid
         
Last update Thur Aug 3 12:20:00 2021
----------------------------------------------------------------------
"""

import numpy as np

def struct2vector(struct, attrname):
    L = len(struct)
    vector = np.empty([L, 1])
    for i, stri in enumerate(struct):
        vector[i] = getattr(struct[i], attrname)

    return vector