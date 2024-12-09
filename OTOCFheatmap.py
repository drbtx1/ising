# -*- coding: utf-8 -*-
"""
Created on Sat Dec  7 10:46:22 2024

@author: dbtx
"""





from Ising_functions import *
from scipy import linalg
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np

J = 1
#h1 = 0.5
#h2 = 0.5
N = 7
h = 1
#g = 0.8
periodic = False
longitudinal_field = sz
transverse_field1 = sx
transverse_field2 = sz
nonHermitianTerms = {0:1j*sy, 6:1j*sy}


def findHamiltonian(N,J,h1,h2,g,periodic,longitudinal_field,transverse_field1,transverse_field2):   
    tensor_sum_J = build_J(N, periodic,longitudinal_field)
    tensor_sum_h1 = build_h(N, periodic, transverse_field1)
    tensor_sum_h2 = build_h(N, periodic, transverse_field2)
    tensor_sum_g = build_g_multiple_tensor(N,periodic,nonHermitianTerms)
    H = -J*tensor_sum_J - h1*tensor_sum_h1 - h2*tensor_sum_h2  +g*tensor_sum_g
    eigenvalues, eigenvectors = np.linalg.eig(H)
    #eigenvectors = np.transpose(eigenvectors)
    return eigenvalues, eigenvectors, H

def V(site,operator):
    #create list of operators
    operators = []
    for index in range(0,N): 
        operators.append(I)
    operators[site] = operator  
    product = operators[0]
    for index in range(1,N):
        product = np.kron(product, operators[index])
    return product 

def W(site,operator):
    #create list of operators
    operators = []
    for index in range(0,N): 
        operators.append(I)
    operators[site] = operator  
    product = operators[0]
    for index in range(1,N):
        product = np.kron(product, operators[index])
    return product 



def U(H, time):
    return linalg.expm(-1j*H*time/h)

#print(eigenvalues)
#t = []
#U = []
#F = []
w = W(1,sz)
v = V(5,sz)
wcc = np.conj(w).T
vcc = np.conj(v).T

#print(w)
#print(v)
#dt = 1
#count = 50
'''for i in range(0,count):
    t.append(i*dt)
    U.append(linalg.expm(-1j*H*t[i]/h))'''
    
#print(t)
#print(U[15])    

def isComplex(x):
    if isinstance(x, complex):
        return True
    else:
        return False
garray = []     
index = 5   
timestep = 0.1
maxtime = 1 
for y in np.arange(0,4,0.01):
    for x in np.arange(0,2,0.01):
        imaginaryTriggered = False
        temparray = []
        for g in np.arange(0,2,0.1):
            eigenvalues, eigenvectors, H = findHamiltonian(N,J,x,y,g,
                                                    periodic,longitudinal_field,
                                                    transverse_field1,transverse_field2)
            psi = eigenvectors[:,index]
            rho = np.outer(psi, np.conj(psi).T)
    
            time = 0 
            
            while imaginaryTriggered == False and time <= maxtime:
                u = U(H, time)
                ucc = np.conj(u).T
                
                F = np.trace((ucc@w@u)@v@rho@(ucc@wcc@u)@vcc)/np.sqrt(np.trace((ucc@w@u)@v@rho@vcc@(ucc@wcc@u))*
                        np.trace(v@(ucc@w@u)@rho@(ucc@wcc@u)@vcc))
                
                if isComplex(F):
                    imaginaryTriggered = True
                    break
                
            if imaginaryTriggered:
                temparray.append(g)
            else: temparray.append(index)    
            
                
    garray.append(temparray)        
                    
    
    print(garray)
                   
sns.heatmap(garray)
#sns.show()    
#plt.colorbar(F)
#plt.plt(F)                
                    
                
   
        
        
        
        #F.append(np.trace((ucc@wcc@u)@vcc@(ucc@w@u)@v@rho))
    #plt.plot(t,F)
    #print(F)
'''F = []
for i in range(0,count):
    wcc = np.conj(w).T
    vcc = np.conj(v).T
    u = U[i]
    ucc = np.conj(u).T
    #F.append(np.trace((ucc@w@u)@vcc@(ucc@wcc@u)@v))
    F.append(np.trace((ucc@w@u)@v@(ucc@w@u)@v))
plt.plot(t,F)'''
plt.show()     
   