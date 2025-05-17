# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 11:54:28 2025

@author: dabuch
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


start = time.time()

#Following lines make plots look a little more Latex-like
#mpl.rcParams['mathtext.fontstyle'] = 'cm' 
mpl.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmss10'
mpl.rcParams['font.monospace'] = 'cmtt10'
#mpl.rcParams["axes.formatter.use_mathtext"] = True # to fix the minus signs

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
                



def findPerturbedHamiltonian(N,J1,J2,h,g, periodic, longitudinal_field1, longitudinal_field2,
                             transverse_field, nonHermitianTerms ):
    tensor_sum_g = build_g_multiple_tensor(N,periodic,nonHermitianTerms)  
    #tensor_sum_g = build_g_single_tensor(N,periodic,site = 0)
    
    tensor_sum_J1 = build_J(N,periodic, longitudinal_field1)
    tensor_sum_J2 = build_J(N,periodic, longitudinal_field2)
    tensor_sum_h = build_h(N,periodic, transverse_field)
    #print(tensor_sum_g)    
    H = -J1*tensor_sum_J1/4 -J2*tensor_sum_J2/4 - h*tensor_sum_h/2 + g*tensor_sum_g 
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


    
J1 = 0.501
J2 = 0.499
h = 0.0
#g = 0.2
#limit = 10e-4

longitudinal_field1 = sx
longitudinal_field2 = sy
transverse_field = sz
N = 4
periodic = False


allUp = []
for n in range(N):
    allUp.append(0)
def makePermutations(oldPerm, site, permutations, N):
    temp = oldPerm.copy()
    temp[site] = 1
    if site == N-1:        
        
        permutations.append(oldPerm)
        permutations.append(temp)
        return
    else:
        makePermutations(oldPerm, site + 1, permutations,N)
        makePermutations(temp, site + 1, permutations,N)
        
def startPermutations(N):      
    
    permutations = []
    makePermutations(allUp,0, permutations,N)   
    return permutations     

def simplifyPermutations(items):
    N = len(items[0])
    reducedList = []
    for item in items:      
        new = True
        for n in range(N):
            temp = item[n:N] + item[0:n]
            #print(temp)
            if temp  in reducedList:
                new = False
        if new == True :   
            reducedList.append(temp)
    #print(reducedList)        
    return reducedList   
        

        
def makeStrings(): 
    strings = []
    permutations = simplifyPermutations(startPermutations(N)) 
    string = ""
    #dictionary hold site of perturbation as key and perturbing operator as value
    for item in permutations: 
        if item != allUp:
            idx = 0
            #string = "nonHermitianTerms = {"
            string = "{"
            for entry in item:
                if entry == 1:
                    string = string + str(idx) + ":splus, "
                idx += 1    
        
            string = string[:-2]
        
            string = string + "}"    
            strings.append(string)
    #print(strings)        
    return strings        

strings = makeStrings()
for string in strings:   
    #print(string)
    nonHermitianTerms = eval(string)
    print(nonHermitianTerms)
    print(type(nonHermitianTerms))
    gamma = []
    realPart = []
    imagPart = []
    for g in np.arange(0,1, 0.01):
        ham = findPerturbedHamiltonian(N,J1,J2,h,g, periodic, 
                                                     longitudinal_field1,longitudinal_field2, transverse_field,
                                                     nonHermitianTerms)
    
    
        ev, vecs = np.linalg.eig(ham)  
        for e in ev:
            realPart.append(re(e))
            gamma.append(g)    #will have g appended multiple times
            if im(e) > 10e-4:
                imagPart.append(im(e))
            else:
                imagPart.append(0)    
            
    #print(ev)    
#realPartList = makeList(realPart, N)  
#imagPartList = makeList(imagPart, N)    
#gammaList = makeList(gamma, N)  
    gamma_array = np.reshape(gamma, (int(len(gamma)/2**N), 2**N)).T
    realPart_array = np.reshape(realPart, (int(len(realPart)/2**N), 2**N)).T
    imagPart_array = np.reshape(imagPart, (int(len(realPart)/2**N), 2**N)).T
# creating an empty canvas
    fig = plt.figure()
 
# defining the axes with the projection
# as 3D so as to plot 3D graphs
    ax = plt.axes(projection="3d")

# plotting a 3D line graph with X-coordinate,
# Y-coordinate and Z-coordinate respectively
    ax.set_xlabel("Gamma")
    ax.set_ylabel("Imaginary part")
    ax.set_zlabel("Real Part")
    ax.set_title("N = " + str(N) + " perturbations: " + string )

    for i in range(0,2**N):
        colors = np.where(imagPart_array[i,:] == 0,  'blue', 'red')
        ax.scatter3D(gamma_array[i,:], imagPart_array[i,:], realPart_array[i,:], marker = ".", c = colors, s = 0.2 )
    path = "C:/Users/dabuch/.spyder-py3/research/Heisenberg/OBC/x501y499for" + str(N) +"spins/"
    
    
    name = ""
    for term in nonHermitianTerms:
        name = string.replace(":","")     
       
            
        
        
    plt.savefig(path +  name + ".png")    
    plt.show()


#print(gamma_array[0])
end = time.time()
print(end-start)