# -*- coding: utf-8 -*-
"""
Generate database from thermo.inp (NASA)

Created on Mon Jun 23 16:08:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import re


def FullName2name(specie):
    name = specie
    if name[-1] == '+':
        name = name[0:-1] + 'plus'
    elif name[-1] == '-':
        name = name[0:-1] + 'minus'
    name = re.sub('[()]', 'b', name)
    name = re.sub('\W', '_', name)
    if re.match('[0-9]', name[0]):
        name = re.sub('[0-9]', 'num_', name[0]) + name

    return name
