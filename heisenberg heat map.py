# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 20:26:17 2025

@author: dbtx
"""




import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import time


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



'''    
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
    return product    '''
                



def findHamiltonian(N,J1,J2,h, periodic, longitudinal_field1, longitudinal_field2,transverse_field ):
    tensor_sum_J1 = build_J(N,periodic, longitudinal_field1)
    tensor_sum_J2 = build_J(N,periodic, longitudinal_field2)
    tensor_sum_h = build_h(N,periodic, transverse_field)
        
    H = -J1*tensor_sum_J1 - J2*tensor_sum_J2 - h*tensor_sum_h
    #vals, vecs = np.linalg.eigh(H)
    
    return H




    
J = 1.00
#g = 0.2
#limit = 10e-4

longitudinal_field1 = sx
longitudinal_field2 = sy
transverse_field = sz
N = 7
periodic = False

Jx = []
Jy = []
minj = -2
maxj = 2
for j in np.arange(minj,maxj,0.1):
    Jx.append(j)
    Jy.append(j)
#print(Jx)
h = 0
ground = np.empty((len(Jx), len(Jy)))
#print(Jx)
for idxjx in range(len(Jx)):
    for idxjy in range(len(Jy)):
        jx = Jx[idxjx]
        jy = Jy[idxjy]
        vals, vecs = np.linalg.eigh(findHamiltonian(N,jx,jy,h, periodic, longitudinal_field1, longitudinal_field2,
                                     transverse_field ))
        ground[idxjx][idxjy] = min(vals)
        


plt.xlabel(r'$J_{x}$')
plt.ylabel(r'$J_{y}$')
plt.title("Ground state of Heisenberg model (h = 0)")
plt.imshow(ground, extent = [minj,maxj,minj,maxj], origin='lower')
plt.colorbar()


#print(gamma_array[0])
end = time.time()
print(end-start)