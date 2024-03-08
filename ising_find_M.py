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


N = 6
periodic = False

site = 1


h_per_J = []
avgExpValueAtSite = []


#build list of operators to find expectation value of nth site
#nth element in each list holds tensor product of sigma operator at nth site
#and all other sites holding identity operator
#x_operator = []
#y_operator = []
operator = []

for site in range(0,N):
    #x_operator.append(expectationOperator(site,N,sx))
    #y_operator.append(expectationOperator(site,N,sy))
    operator.append(expectationOperator(site,N,sz))
 
maxh = 100   
for h in range(0,maxh): 
    #hold expectation values at nth site 
    exp = []    
    eigenvalues, eigenvectors = findHamiltonian(N,J,h/maxh,periodic)
    exp = findExpectationValues(N, site,operator, eigenvectors)
    avg = sum(exp)/N   
    h_per_J.append((h/maxh)/J)
    avgExpValueAtSite.append(avg)
    
plt.plot(h_per_J, avgExpValueAtSite)
plt.show()    



    
   