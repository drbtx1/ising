# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 12:24:35 2024

@author: dabuch
"""

#Code to find Hamiltonian and calculate expectation valuesfor 1-d Ising model.
# z interactions occur between nearest neighbors, and sites have individual 
#x interactions. Periodic boundary conditions can be turned on or off. Code 
#assumes last site can have z interaction 
#without periodic boundary conditions
from Ising_functions import *
import numpy as np
from matplotlib import pyplot as plt
#import csv
import pandas as pd

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])


J = 1.0
maxh = 4
steps = 10000
#N = 3
periodic = True

#site = 1
#eV = 0  #eigenvector to find M for (average sz across all sites in chain)







for N in range(3,8):
    for eV in range(0, 2**N):
        findAndSaveMagnetization(N,periodic,eV,J,maxh,steps)

   





