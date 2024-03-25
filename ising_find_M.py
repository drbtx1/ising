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

longitudinal_field = sz
transverse_field = sx
direction_of_M = sz


J = 1.0
maxh = 10
steps = 1000
#N = 3
periodic = False

#site = 1
#eV = 0  #eigenvector to find M for (average sz across all sites in chain)







for N in range(4,5):
    for eV in range(0, 1): #2**N):
        findAndSaveMagnetization(N,periodic,eV,J,maxh,steps,longitudinal_field,transverse_field, direction_of_M )

   





