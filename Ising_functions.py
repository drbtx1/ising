# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:12:29 2024

@author: dabuch
"""

import numpy as np

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
splus = (sx + 1j*sy)/2
sminus = (sx - 1j*sy)/2


#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing nearest neighbor interactions
def build_J(N,periodic):
    tensor_sum_J = 0
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
    return tensor_sum_J        
            
#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing localized term in x directions at each site
def build_h(N, periodic):
    #intra-site   
    tensor_sum_h = 0
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
    return tensor_sum_h 

def build_g(N,periodic,nonHermitianTerms):
    #perturbed sites 
    #create list of operators, set all to identity, change terms that should be
    #non-Hermitian perturbations, then find kronecker product 
    operators = []
    for index in range(0,N): 
        operators.append(I)
    for site in nonHermitianTerms:
        operators[site] = nonHermitianTerms[site]
    #print(operators)    
    product = operators[0]
    for index in range(1,N):
        product = np.kron(product, operators[index])
    tensor_sum_g = product 
    #print(product)
    return tensor_sum_g    
    #print(nonHermitianTerms)
    
            
#takes number of sites N, energy scaling J and h, and boolean periodic (boundary conditions)
#return lists of eigenvalues and eigenfunctions for corresponding Hamiltonian 
def findHamiltonian(N,J,h, periodic):
   

    H = J*tensor_sum_J + h*tensor_sum_h  
    #print(H)

    eigenvalues, eigenvectors = np.linalg.eig(H)
    #print(eigenvalues)        
    return eigenvalues, eigenvectors



def findPerturbedHamiltonian(N,J,h,g, periodic, nonHermitianTerms ):
    tensor_sum_g = build_g(N,periodic,nonHermitianTerms)    
    tensor_sum_J = build_J(N,periodic)
    tensor_sum_h = build_h(N,periodic)
    
    H = J*tensor_sum_J + h*tensor_sum_h + g*tensor_sum_g 
    

    eigenvalues, eigenvectors = np.linalg.eig(H)
    #print(g*tensor_sum_g)        
    return eigenvalues, eigenvectors




#takes N total sites, and nth specific site, and a Pauli operator
#returns tensor product of pauli on n and identity operator on all other sites
#this is used to expectation value of spin in various directions on site n 
def expectationOperator(n,N, pauli):      #N sites, nth location
    operator = []
    #opz = []
    for i in range(N):
        operator.append(I)
        #opz.append(I)
    operator[n] = pauli
    #opz[n] = sz    
    product = operator[0]
    #zproduct = opz[0]
    for index in range(1,N):
        product = np.kron(product, operator[index])
        #zproduct = np.kron(zproduct, opz[index])
    return product   #, zproduct

#expects a list of tensor products (give expectation values at each site), and a list of eigenvectors  
#returns a list of expectation values (corresponding to each eigenvector) for that site
def findExpectationValues(N,site,pauliOperatorProduct, eigenvectors):
    
    ExpValAtSite = []
    #print("site = " + str(site))
           
    for v in range(0,len(eigenvectors)):           
        x = np.matmul(np.matmul(eigenvectors[v].conj().T,pauliOperatorProduct[site]), eigenvectors[v])
        #print(x)
        ExpValAtSite.append(x)
            
    return ExpValAtSite        
            
                    
        