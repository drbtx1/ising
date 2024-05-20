# -*- coding: utf-8 -*-
"""
Created on Sat May 11 15:13:01 2024

@author: dbtx
"""

from Ising_functions import *
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
#import csv
import pandas as pd
from sympy import *

h = symbols("h")




I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])

longitudinal_field = sz
transverse_field = sx
direction_of_M = sx


J = 1
maxh = 0.0000000000000000000001
steps = 2
N = 3
periodic = False
h = 0.00001

H = -J*build_J(N,periodic, longitudinal_field) - h*build_h(N, periodic, transverse_field)
print(H)
#print(Matrix(H).eigenvects())
'''lamda = symbols('lamda')
#p = Array(H).charpoly(lamda)
#factor(p.as_expr())'''
'''for step in range(steps):
    h = step*maxh/steps
    print("h = " + str(h))
    eigenvalues, eigenvectors = findHamiltonian(N,J,h, periodic,longitudinal_field,transverse_field)
    idx = eigenvalues.argsort()[::-1]   
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:,idx]
    print(eigenvectors)
    print(eigenvalues)'''
    
eigenvalues, eigenvectors = findHamiltonian(N,J,h, periodic,longitudinal_field,transverse_field)
idx = eigenvalues.argsort()[::-1]   
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:,idx]
print(eigenvalues)
print(eigenvectors)
#print([list(tup[2][0]) for tup in H.eigenvects()])
