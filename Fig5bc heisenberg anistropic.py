# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 13:40:13 2025

@author: dbtx
"""


import numpy as np
#import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
from sympy import Matrix, re, im
#import seaborn as sns 

start = time.time()

#Following lines make plots look a little more Latex-like
#mpl.rcParams['mathtext.fontstyle'] = 'cm' 
mpl.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmss10'
mpl.rcParams['font.monospace'] = 'cmtt10'
mpl.rcParams["axes.formatter.use_mathtext"] = True # to fix the minus signs

#I = Matrix([1,0],[0,1])
I = np.eye(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
splus = (sx + 1j*sy)/2
sminus = (sx - 1j*sy)/2


#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing nearest neighbor interactions
def build_J(N,periodic, longitudinal_field):
    J_operator = longitudinal_field
    #tensor_sum_J = 0
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
            if dim == 0:
                tensor_sum_J = product
            else:
                tensor_sum_J = tensor_sum_J + product
    return tensor_sum_J        
            
#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing localized term in  transverse direction at each site
def build_h(N, periodic, transverse_field):
    #intra-site   
    h_operator = transverse_field
    #tensor_sum_h = 0
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)
        operators[dim] = h_operator  
        product = operators[0]
        for index in range(1,N):
            product = np.kron(product, operators[index])
        if dim == 0:
            tensor_sum_h = product
        else:    
            tensor_sum_h = tensor_sum_h + product    
    return tensor_sum_h 



    
#Takes a dictionary of sites and perturbations, call find_perturbation to 
#find sum of kronecker product of each perturbed site 
def build_g_multiple_tensor(N,periodic,nonHermitianTerms):  
    tensor_sum_g = np.zeros((2**N, 2**N))
    for site in nonHermitianTerms:
        #if site == 0:
        #    tensor_sum_g = find_perturbation(N,site,nonHermitianTerms[site]) 
        #else:
        tensor_sum_g = tensor_sum_g + find_perturbation(N,site,nonHermitianTerms[site])         
    return tensor_sum_g  

def build_g_samesite(N,periodic,site):  
    tensor_sum_g = np.zeros((2**N, 2**N))
    for site in nonHermitianTerms:
        #if site == 0:
        #    tensor_sum_g = find_perturbation(N,site,nonHermitianTerms[site]) 
        #else:
        tensor_sum_g = tensor_sum_g + find_perturbation(N,site,nonHermitianTerms[site])         
    return tensor_sum_g     
    
def build_g_single_tensor(N,periodic,site):   #handle p = q
    #perturbed sites 
    #create list of operators, set all to identity, change terms that should be
    #non-Hermitian perturbations, then find kronecker product 
    operators = []
    for index in range(0,N): 
        operators.append(I)
    operators[site] = sx
    
    #print(operators)    
    product = operators[0]
    for index in range(1,N):
        product = np.kron(product, operators[index])
    tensor_sum_g = product 
    #print(product)
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
                



'''def findPerturbedHamiltonian(N,J,h,g, periodic, longitudinal_field,
                             transverse_field, nonHermitianTerms, ):
    tensor_sum_g = build_g_multiple_tensor(N,periodic,nonHermitianTerms)  
    #tensor_sum_g = build_g_single_tensor(N,periodic,site = 0)
    
    tensor_sum_J = build_J(N,periodic, longitudinal_field)
    tensor_sum_h = build_h(N,periodic, transverse_field)
    #print(tensor_sum_g)    
    H = -J*tensor_sum_J/4 - h*tensor_sum_h/2 + g*tensor_sum_g 
    return(H)'''

def findPerturbedHamiltonian(N,J,deltaJ,h,g, periodic, longitudinal_field1, longitudinal_field2,
                             transverse_field, nonHermitianTerms ):
    tensor_sum_g = build_g_multiple_tensor(N,periodic,nonHermitianTerms)  
    #tensor_sum_g = build_g_single_tensor(N,periodic,site = 0)
    
    tensor_sum_J1 = build_J(N,periodic, longitudinal_field1) + build_J(N,periodic, longitudinal_field2)
    tensor_sum_J2 = build_J(N,periodic, longitudinal_field1) - build_J(N,periodic, longitudinal_field2)
    tensor_sum_h = build_h(N,periodic, transverse_field)
    #print(tensor_sum_g)    
    H = -J*tensor_sum_J1/4 -deltaJ*tensor_sum_J2/4 - h*tensor_sum_h/2 + g*tensor_sum_g 
    return(H)




    

gList = []
for gamma in np.arange(0,2, 0.01):
    gList.append(gamma)

J = 1    
deltaJ = 0
h = 0.1
#g = 0.2
#limit = 10e-4

longitudinal_field1 = sx
longitudinal_field2 = sy
transverse_field = sz
N = 4
periodic = False

#dictionary hold site of perturbation as key and perturbing operator as value
#nonHermitianTerms = {0:splus + sminus}






p = 0
q = 1
nonHermitianTerms = {p: splus, q: -1*sminus}
realPart = []
imagPart = []
purelyReal = []
for idx in range(len(gList)):
    g = gList[idx]
    ham = findPerturbedHamiltonian(N,J,deltaJ,h,g, periodic, 
                                             longitudinal_field1, longitudinal_field2, transverse_field,
                                             nonHermitianTerms)
    ev, vecs = np.linalg.eig(ham)
    tempReal = []
    tempImag = []
    tempPurelyReal = []
    for e in ev:
        if abs(np.imag(e)) < 10e-5:
            tempImag.append(np.nan)
            tempReal.append(np.nan) 
            tempPurelyReal.append(np.real(e))
        
        
        else:           
            tempImag.append(np.imag(e))
            tempReal.append(np.real(e))
            tempPurelyReal.append(np.nan)
    realPart.append(tempReal)
    imagPart.append(tempImag)
    purelyReal.append(tempPurelyReal)
        

realPartArray = np.array(realPart)
imagPartArray = np.array(imagPart)
purelyRealArray = np.array(purelyReal) 
    
for n in range(0, 2**N):
    #plt.plot(gList, realPartArray[:,n])
    #plt.scatter(gList, purelyRealArray[:,n],marker = ".", color = "blue", s=0.03 )    
    #plt.scatter(gList, realPartArray[:,n],marker = ".", color = "red", s=0.03 )
    plt.scatter(gList, imagPartArray[:,n],marker = ".", color = "red", s=0.03 )
#print(gList)    
            
  
        
            
            
            
                    
                        
                
    
    

#plt.title("PT Threshold")
#plt.legend()
#plt.xscale("log")
plt.xlabel(r'$\gamma_{PT}/J/$')
plt.ylabel(r'$Re(E)/J$')
#plt.yticks([1,1.5,2])
#plt.xticks([10**-1,10**0,10**1])'''
#plt.ylim(0,0.15)
#plt.xlim([min,max])
#plt.ylim([min,max])                                       




plt.show()


#print(gamma_array[0])
end = time.time()
print(end-start)