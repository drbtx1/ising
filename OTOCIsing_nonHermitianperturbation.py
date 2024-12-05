# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 19:48:25 2024

@author: dbtx
"""



from Ising_functions import *
from scipy import linalg

J = 1
h1 = 0.1
h2 = 0.1
N = 6
h = 1
g = 0.01
periodic = False
longitudinal_field = sz
transverse_field1 = sx
transverse_field2 = sz
nonHermitianTerms = {0:1j*sy}


def findHamiltonian(N,J,h1,h2,periodic,longitudinal_field,transverse_field1,transverse_field2):   
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

eigenvalues, eigenvectors, H = findHamiltonian(N,J,h1,h2,
                                            periodic,longitudinal_field,
                                            transverse_field1,transverse_field2)

#print(eigenvalues)
t = []
U = []
F = []
w = W(0,sx)
v = V(5,sz)
print(w)
print(v)
dt = 1
count = 500
for i in range(0,count):
    t.append(i*dt)
    U.append(linalg.expm(-1j*H*t[i]/h))
    
#print(t)
#print(U[15])    

for index in range(0,1):
    psi = eigenvectors[:,index]
    rho = np.outer(psi, np.conj(psi).T)
    F = []
    for i in range(0,count):
        wcc = np.conj(w).T
        vcc = np.conj(v).T
        u = U[i]
        ucc = np.conj(u).T
        F.append(np.trace((ucc@wcc@u)@vcc@(ucc@w@u)@v@rho))
    plt.plot(t,F)
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
   