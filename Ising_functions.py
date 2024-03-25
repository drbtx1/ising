# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:12:29 2024

@author: dabuch
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
splus = (sx + 1j*sy)/2
sminus = (sx - 1j*sy)/2


#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing nearest neighbor interactions
def build_J(N,periodic, longitudinal_field):
    J_operator = longitudinal_field
    tensor_sum_J = 0
    #cross-site interactions
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)  #set all to identity and adjust appropriates sits below
     #for each iteration of loop, except last one, set that site and the next one to J_operator         
        if dim != N-1:           
            operators[dim] = J_operator
            operators[dim+1] = J_operator
    #only set last site and the next (first) one to  if the periodic bc are needed        
        if dim == N-1 and periodic == True:
            operators[dim] = J_operator
            operators[0] = J_operator
            
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
#describing localized term in  transverse direction at each site
def build_h(N, periodic, transverse_field):
    #intra-site   
    h_operator = transverse_field
    tensor_sum_h = 0
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)
        operators[dim] = h_operator  
        product = operators[0]
        for index in range(1,N):
            product = np.kron(product, operators[index])
        tensor_sum_h = tensor_sum_h + product    
    return tensor_sum_h 


def build_g_single_tensor(N,periodic,nonHermitianTerms):
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
    
    
#Takes a dictionary of sites and perturbations, call find_perturbation to 
#find sum of kronecker product of each perturbed site 
def build_g_multiple_tensor(N,periodic,nonHermitianTerms):  
    tensor_sum_g = 0
    
    for site in nonHermitianTerms:
        tensor_sum_g = find_perturbation(N,site,nonHermitianTerms[site])         
    return tensor_sum_g    
    
    
def find_perturbation(N,site, nonHermitianTerm):
#create list of operators, set all to identity, change term that should be
#non-Hermitian perturbation, then return kronecker product     
    operators = []
    for index in range(0,N): 
        operators.append(I)
    operators[site] = nonHermitianTerm
    product = operators[0]
    for index in range(1,N):
        product = np.kron(product, operators[index]) 
    #print(product)
    return product    
                
#takes number of sites N, energy scaling J and h, and boolean periodic (boundary conditions)
#return lists of eigenvalues and eigenfunctions for corresponding Hamiltonian 
def findHamiltonian(N,J,h, periodic,longitudinal_field,transverse_field):   
    tensor_sum_J = build_J(N, periodic,longitudinal_field)
    tensor_sum_h = build_h(N, periodic, transverse_field)
    H = J*tensor_sum_J + h*tensor_sum_h  
    eigenvalues, eigenvectors = np.linalg.eig(H)
    return eigenvalues, eigenvectors



def findPerturbedHamiltonian(N,J,h,g, periodic, nonHermitianTerms ):
    tensor_sum_g = build_g_single_tensor(N,periodic,nonHermitianTerms)    
    tensor_sum_J = build_J(N,periodic)
    tensor_sum_h = build_h(N,periodic)    
    H = J*tensor_sum_J + h*tensor_sum_h + g*tensor_sum_g     
    eigenvalues, eigenvectors = np.linalg.eig(H)
    return eigenvalues, eigenvectors




#takes N total sites, and nth specific site, and a Pauli operator
#returns tensor product of pauli on n and identity operator on all other sites
#this is used to expectation value of spin in various directions on site n 
def expectationOperator(n,N, pauli):      #N sites, nth location
    operator = []
    for i in range(N):
        operator.append(I)
        
    operator[n] = pauli
        
    product = operator[0]
    
    for index in range(1,N):
        product = np.kron(product, operator[index])
        
    return product   

#expects a list of tensor products (give expectation values at each site), and a list of eigenvectors  
#returns a list of expectation values (corresponding to each eigenvector) for that site
def findExpectationValue(pauliOperatorProduct, eigenvector):
    return np.matmul(np.matmul(eigenvector.conj().T,pauliOperatorProduct), eigenvector)     


def findAndSaveMagnetization(N,periodic,eV,J,maxh,steps,longitudinal_field, transverse_field, direction_of_M):
    h_per_J = []

    #build list of operators to find expectation value of nth site
    #nth element in each list holds tensor product of sigma operator at nth site
    #and all other sites holding identity operator
    #x_operator = []
    #y_operator = []
    operator_list = []
    avgExpValue = []
    realPart = []
    imagPart = []


    for site in range(0,N):
        operator_list.append(expectationOperator(site,N,direction_of_M))
     
       
    for step in range(steps): 
        #hold expectation values at nth site 
        expectation_at_site = [] 
        h = maxh*step/steps
        eigenvalues, eigenvectors = findHamiltonian(N,J,h,periodic,longitudinal_field, transverse_field)
        for site in range(N):
            expectation_at_site.append(findExpectationValue(operator_list[site], eigenvectors[eV]))
        avg = sum(expectation_at_site)/N   
        #print(avg)
        h_per_J.append((h)/J)
        avgExpValue.append(avg)  #average across N sites, all with same eigenvector
        realPart.append(avg.real)
        imagPart.append(avg.imag)
        
        
        
    #data = avgExpValue    
    title = "N = " +str(N)+" , eigenvector = " + str(eV) + ", periodic = " + str(periodic)

    '''filetype = '.csv'
    dest = 'C:/Users/dabuch/Ising Data/Magnetization/Periodic/N' + str(N) + '/' + title + filetype

    data = {"h/J": h_per_J, "Expectation Value averaged over all sites": avgExpValue}
    df = pd.DataFrame(data)
    df.to_csv(dest)'''
    
    #plt.xlabel("h/J")
    plt.ylabel("M (spin averaged over all sites)")
    plt.title(title)   
    #plt.legend(["Real","Imaginary"])
    plt.plot(h_per_J, realPart)
    #plt.plot(h_per_J, imagPart)
    '''image_filetype = '.png'
    dest = 'C:/Users/dabuch/Ising Data/Magnetization/Periodic/N' + str(N) + '/' + title + image_filetype
    plt.savefig(dest)
    plt.show()'''

    
            
                    
        