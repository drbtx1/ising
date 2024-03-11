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

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])


J = 1.0
maxh = 0.1
steps = 100
N = 6
periodic = False

#site = 1
eV = 1  #eigenvector to find M for (average sz across all sites in chain)

h_per_J = []

#build list of operators to find expectation value of nth site
#nth element in each list holds tensor product of sigma operator at nth site
#and all other sites holding identity operator
#x_operator = []
#y_operator = []
operator_list = []
avgExpValue = []


for site in range(0,N):
    #x_operator.append(expectationOperator(site,N,sx))
    #y_operator.append(expectationOperator(site,N,sy))
    operator_list.append(expectationOperator(site,N,sz))
 
   
for step in range(steps): 
    #hold expectation values at nth site 
    expectation_at_site = [] 
    h = maxh*step/steps
    eigenvalues, eigenvectors = findHamiltonian(N,J,h,periodic)
    for site in range(N):
        expectation_at_site.append(findExpectationValue(operator_list[site], eigenvectors[eV]))
    avg = sum(expectation_at_site)/N   
    #print(avg)
    h_per_J.append((h)/J)
    avgExpValue.append(avg)  #average across N sites, all with same eigenvector
title = "N = " +str(N)+" , eigenvector = " + str(eV) + ", periodic = " + str(periodic)
plt.xlabel("h/J")
plt.ylabel("M (spin averaged over all sites)")
plt.title(title)     
plt.plot(h_per_J, avgExpValue)
plt.show()

   





