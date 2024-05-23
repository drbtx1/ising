# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:51:49 2024

@author: dbtx
"""

#Code to find Hamiltonian and calculate expectation valuesfor 1-d Ising model.
#z interactions occur between nearest neighbors of strength J, and sites have 
#individual x interactions with applied magnetic field, h.
#Periodic boundary conditions can be turned on or off. 
from Ising_functions import *
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
#import csv
import pandas as pd

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])

longitudinal_field = sz
transverse_field = sx
direction_of_M = sx


J = 1
maxh = 3
steps = 200
N = 4
periodic = False

h_per_J = []

for step in range(steps):
    h = maxh*step/steps
    h_per_J.append(h/J)

def geteigenvalues_vs_h(N,periodic,eV,J,maxh,steps,longitudinal_field, transverse_field):
    eV_list = []
    for step in range(steps):
        h = maxh*step/steps
        eigenvalues, eigenvectors = findHamiltonian(N,J,h,periodic,longitudinal_field, transverse_field)
        #assure sorted in decreasing order
        idx = eigenvalues.argsort()[::-1]   
        eigenvalues = eigenvalues[idx]
        eV_list.append(eigenvalues[eV])
    plt.plot(h_per_J, eV_list)
    return eV_list    

#print(geteigenvalues_vs_h(3,periodic,0,J,maxh,steps,longitudinal_field, transverse_field))
for vec in range(0,2**N):
    geteigenvalues_vs_h(N,periodic,vec,J,maxh,steps,longitudinal_field, transverse_field)
#mpl.rcParams['text.usetex'] = False
'''mpl.use('Agg')
dest = 'C:/Users/dbtx/.spyder-py3/eigs.png'

plt.title("Plot of eigenvalue vs transverse field, " + str(N) + " spins")  '''
#plt.savefig(dest) 
plt.xlabel("h/J")
plt.ylabel("Energy") 
plt.show()    

        
        