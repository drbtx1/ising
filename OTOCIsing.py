# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 11:03:32 2024

@author: dbtx
"""

from Ising_functions import *
from scipy import linalg

J = 1
h1 = 0.1
h2 = 0.1
N = 3
h = 1
periodic = False
longitudinal_field = sz
transverse_field1 = sx
transverse_field2 = sz


def findHamiltonian(N,J,h1,h2,periodic,longitudinal_field,transverse_field1,transverse_field2):   
    tensor_sum_J = build_J(N, periodic,longitudinal_field)
    tensor_sum_h1 = build_h(N, periodic, transverse_field1)
    tensor_sum_h2 = build_h(N, periodic, transverse_field2)
    H = -J*tensor_sum_J - h1*tensor_sum_h1 - h2*tensor_sum_h2 
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

eigenvalues, eigenvectors, H = findHamiltonian(N,J,h1,h2,
                                            periodic,longitudinal_field,
                                            transverse_field1,transverse_field2)

#print(eigenvalues)
t = []
U = []
F = []
w = W(0,sx)
v = V(2,sz)
dt = 0.1
count = 500
for i in range(0,count):
    t.append(i*dt)
    U.append(linalg.expm(-1j*H*t[i]/h))
    
#print(t)
#print(U[15])    

'''for index in range(0,2**N):
    psi = eigenvectors[:,index]
    rho = np.outer(psi, np.conj(psi).T)
    F = []
    for i in range(0,count):
        wcc = np.conj(w).T
        vcc = np.conj(v).T
        u = U[i]
        ucc = np.conj(u).T
        F.append(np.trace((ucc@w@u)@vcc@(ucc@wcc@u)@v@rho))
    plt.plot(t,F)
    print(F)'''
F = []
for i in range(0,count):
    wcc = np.conj(w).T
    vcc = np.conj(v).T
    u = U[i]
    ucc = np.conj(u).T
    F.append(np.trace((ucc@w@u)@vcc@(ucc@wcc@u)@v))
plt.plot(t,F)
plt.show()     
   