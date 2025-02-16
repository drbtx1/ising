# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 14:07:41 2025

@author: dbtx
"""



import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
#from mpl_toolkits import mplot3d
#from sortev import *
#from sortcolumns import *
#from sortvecs import *
#from scipy.linalg import eig
from sympy import Matrix, re, im
#from sympy.physics.quantum import TensorProduct
#from sympy import eye
#from sympy import zeros
#import multiprocessing as mp
import seaborn as sns 

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
                



def findPerturbedHamiltonian(N,J,h,g, periodic, longitudinal_field,
                             transverse_field, nonHermitianTerms, ):
    tensor_sum_g = build_g_multiple_tensor(N,periodic,nonHermitianTerms)  
    #tensor_sum_g = build_g_single_tensor(N,periodic,site = 0)
    
    tensor_sum_J = build_J(N,periodic, longitudinal_field)
    tensor_sum_h = build_h(N,periodic, transverse_field)
    #print(tensor_sum_g)    
    H = -J*tensor_sum_J/4 - h*tensor_sum_h/2 + g*tensor_sum_g 
    return(H)

def makeList(M, N):
    list = []
    columns = 2**N
    
    rows = int(len(M)/columns)
    #print(columns)
    #print(rows)
    #print(len(M))
    for row in range(rows):
        temp = []
        for column in range(columns):
            #print(row)
            temp.append(M[row*columns + column])
        list.append(temp)    
    return list


    
J = 1.0
h = 0.0
#g = 0.2
#limit = 10e-4

longitudinal_field = sx
transverse_field = sz
N = 6
periodic = False

#dictionary hold site of perturbation as key and perturbing operator as value
#nonHermitianTerms = {0:splus + sminus}

#gamma = []
#realPart = []
#imagPart = []


p = 0
dg = 0.01
max = 1 
min = -1 
maxIm = np.zeros((int((max-min)/dg),int((max-min)/dg)))
print(maxIm)

for gplus in np.arange(min,max+dg ,dg):
    print(gplus)
    for gminus in np.arange(min,max+dg , dg):
        print(gminus)
        nonHermitianTerms = {p: (gplus + gminus)*sx + 1j*(gplus - gminus)*sy}
        ham = findPerturbedHamiltonian(N,J,h,1, periodic, 
                                                     longitudinal_field, transverse_field,
                                                     nonHermitianTerms)
        ev, vecs = np.linalg.eig(ham)
        #print(ev)
        for e in ev:
            if im(e) > maxIm[int((gplus-min)/dg)][int((gminus-min)/dg)]:
                maxIm[int((gplus-min)/dg)][int((gminus-min)/dg)] = im(e)
            
                    
                
    
    
plt.imshow(maxIm,extent = [-1,1,-1,1], origin='lower')
plt.title("Maximum imaginary part, p = " + str(p))
#plt.xlim([min,max])
#plt.ylim([min,max])
'''ax = sns.heatmap(maxIm)
ax.invert_yaxis()
ax.set_yticks(np.arange(min,max,0.1)) 
#ax1.set_xticks(range(9, 200, 10))
#ax.set_xticklabels(f'{c:.1f}' for c in np.arange(-1.00, 1.00, 0.1))
#ax.set_yticklabels(f'{c:.1f}' for c in np.arange(-1.00, 1.00, 0.1))'''
plt.show()


#print(gamma_array[0])
end = time.time()
print(end-start)