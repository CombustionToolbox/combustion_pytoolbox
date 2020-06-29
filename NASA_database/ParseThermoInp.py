# -*- coding: utf-8 -*-
"""
Generate database from thermo.inp (NASA)

Created on Mon Jun 23 16:08:00 2020

@author: Alberto Cuadra Lara
         PhD Candidate - Group Fluid Mechanics
         Office 1.1.D17, Universidad Carlos III de Madrid
"""
import re
import pickle
import numpy as np
from NASA_database.FullName2name import FullName2name

def ParseThermoInp(reducedDB):
    filepath = 'Databases/strMaster_reduced.pkl'
    if not existFile(filepath):
        fid = open('D:/Google Drive/Phd/Combustion_Toolbox/ThermochemicalCode_Python/ThermochemicalCode_Python/NASA_database/thermo.inp', 'r')
        
        if reducedDB:
            print('Loading Reduced NASA database ...', end = " ")
        else:
            print('Loading NASA database ...', end = " ")
        struct = {}
        while True: 
            tline = fid.readline()
            if tline == '':
                break
            if tline[0] == '!':
                continue
            if 'thermo' in tline:
                tline = fid.readline()
                continue
            if 'END' in tline:
                continue
            
            aux = StrMaster(tline, fid)
            struct.update({aux.name: aux})

        if reducedDB:
            struct = StrMaster_reduced(struct)
            # Save StrMaster reduced 
            f = open("Databases/strMaster_reduced.pkl","wb")
            pickle.dump(struct,f)
            f.close()
            
        print('OK!')
            
        
    else: # Load StrMaster reduced
        with open(filepath, 'rb') as f:
            struct = pickle.load(f)
        print('NASA reduced database loaded from %s' % filepath)
    
    return struct

def StrMaster_reduced(self):
    NameSpecies = list(self)
    NSpecies = len(NameSpecies)
    ind = np.ones((NSpecies, 1))
    patterns = ['plus','minus','AL','Ag','F','CL','B','Ca','I','K',\
    'Li','M','D','S','Rb','Pb','V','W','Z','G','T','Cd','Co','Cr',\
    'Cs','Cu','Ni','U','Na','Nb','Hg','CP','HP']
    
    for key in NameSpecies:
        if key.startswith('P'):
            del self[key]
            continue
        for pattern in patterns:
            if pattern in key:
                del self[key]
                break
            
    return self

def existFile(filepath):
    import os 
    return os.path.exists(filepath)

class StrMaster():
    def __init__(self, tline, fid):
        self.struct(tline, fid)
    def struct(self, tline, fid):
        self.FullName = tline[0:16].replace(' ', '')
        self.name = FullName2name(self.FullName)
        self.comments = tline[18::]
        tline = fid.readline()
        self.ctTInt = int(tline[0:2])
        self.txRefCode = tline[3:9]
        self.txFormula = tline[10:50]
        self.swtCondensed = int(tline[50:52])
        self.mm = float(tline[52:65])
        self.Hf0 = float(tline[65:80])
        
        if not self.ctTInt:
            tline = fid.readline()
            self.tRange = [float(i) for i in tline[0:22].split()]
            self.tExponents = [float(i) for i in tline[23:63].split()]
            self.Hf298De10 = [float(i) for i in tline[65::].split()]
        else:
            self.tRange = []
            self.tExponents = []
            self.Hf298De10 = []
            self.a = []
            self.b = []
            for i in range(0, self.ctTInt):
                tline = fid.readline()
                self.tRange.append([float(i) for i in tline[0:22].split()]) 
                self.tExponents.append([float(i) for i in tline[23:63].split()])
                self.Hf298De10.append([float(i) for i in tline[65::].split()])
                tline = fid.readline()
                a1 = float(tline[0:16])
                a2 = float(tline[16:32])
                a3 = float(tline[32:48])
                a4 = float(tline[48:64])
                a5 = float(tline[64:80])
                tline = fid.readline()
                a6 = float(tline[0:16])
                a7 = float(tline[16:32])
                a8 = 0.
                b1 = float(tline[48:64])
                b2 = float(tline[64:80])
                self.a.append([a1, a2, a3, a4, a5, a6, a7, a8])
                self.b.append([b1, b2])
