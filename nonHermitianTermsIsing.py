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
import scipy as sp
from matplotlib import pyplot as plt

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])


J = 1.0
h = 0.0
g = 0.5
limit = 10e-4

longitudinal_field = sx
transverse_field = sz
N = 7
periodic = False

#dictionary hold site of perturbation as key and perturbing operator as value
nonHermitianTerms = {2:splus, 5 :sminus}
#for k, v in nonHermitianTerms.items():
#    print(type(k))
#    print(type(v))
#eigenvalues, eigenvectors = findPerturbedHamiltonian(N,J,h,g, periodic, 
                                                     #longitudinal_field, transverse_field,
                                                     #nonHermitianTerms)
    
print(findPerturbedHamiltonian(N,J,h,g, periodic, 
                                                     longitudinal_field, transverse_field,
                                                     nonHermitianTerms))    

#eigenvalues, eigenvectors = findHamiltonian(N,J,h, periodic,longitudinal_field,
#                                            transverse_field)    
#print(eigenvectors)

'''for element in eigenvalues:
    if np.imag(element) < limit:
        element = np.real(element)
    print(element)    ''' 
print(eigenvalues)

