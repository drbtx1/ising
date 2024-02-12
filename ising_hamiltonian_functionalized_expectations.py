# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 14:37:39 2024

@author: dabuch
"""

# -*- coding: utf-8 -*-
#Code to find Hamiltonian for 1-d Ising model. z interactions occur between nearest neighbors, and sites have individual 
#x interactions. Periodic boundary conditions can be turned on or off. Code assumes last site can have z interaction 
#without periodic boundary conditions

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import timeit

periodic = False

#N = 4   #specify number of sites
#code does not take into account special case of N=2 for periodic boundary conditions

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])

J = 1.0
h = 0.01

N = 9

def findHamiltonian(N):
    tensor_sum_J = 0
    tensor_sum_h = 0
    #total_tensor = 0
    
    #cross-site interactions
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)  #set all to identity and adjust appropriates sits below
     #for each iteration of loop, except last one, set that site and the next one to sz         
        if dim != N-1:           
            operators[dim] = sz
            operators[dim+1] = sz
    #only set last site and the next (first) one to sz if the periodic bc are needed        
        if dim == N-1 and periodic == True:
            operators[dim] = sz
            operators[0] = sz
            
        #create tensor product of operators     
        product = operators[0]
        for index in range(1,N):
            product = np.kron(product, operators[index])
    #check not to add extraneous identity in the case that loop is in last iteration
    #and periodic bc are not needed        
        if dim != N-1 or periodic == True:     
            tensor_sum_J = tensor_sum_J + product


    #intra-site   
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)
        operators[dim] = sx  
        product = operators[0]
        for index in range(1,N):
            product = np.kron(product, operators[index])
        tensor_sum_h = tensor_sum_h + product    
    #print(tensor_sum_h)  

    H = J*tensor_sum_J + h*tensor_sum_h  
    #print(H)

    eigenvalues, eigenvectors = np.linalg.eig(H)
    #print(eigenvalues)        
    return eigenvalues, eigenvectors

eigenvalues, eigenvectors = findHamiltonian(N)
#print(eigenvalues)

exp_x = []
exp_z = []
vec = []

def buildOps(n,N):      #N sites, nth location
    opx = []
    opz = []
    for i in range(N):
        opx.append(I)
        opz.append(I)
    opx[n] = sx
    opz[n] = sz    
    xproduct = opx[0]
    zproduct = opz[0]
    for index in range(1,N):
        xproduct = np.kron(xproduct, opx[index])
        zproduct = np.kron(zproduct, opz[index])
    return xproduct, zproduct   

for loc in range(0,N):
    vec.append(N)
    #buildup tensor product as function
    op_x, op_z = buildOps(loc,N)    
    
    
    x = np.matmul(np.matmul(eigenvalues[loc],op_x), eigenvalues[loc])
    exp_x.append(x)
    
    z = np.matmul(np.matmul(eigenvalues[loc],op_z), eigenvalues[loc])
    exp_z.append(y)
    
    
    

              
    


