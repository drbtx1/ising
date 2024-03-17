# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 15:31:16 2024

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

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])


J = 1
h = 0.0
g = 0.6

N = 5
periodic = False

#dictionary hold site of perturbation as key and perturbing operator as value
nonHermitianTerms = {0:splus, 2:sminus}

eigenvalues, eigenvectors = findPerturbedHamiltonian(N,J,h,g, periodic, nonHermitianTerms)

print(eigenvectors)
print(eigenvalues)

