# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:18:58 2025

@author: dbtx
"""



import numpy as np
#import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
#from mpl_toolkits import mplot3d
#from sortev import *
#from sortcolumns import *
#from sortvecs import *
#from scipy.linalg import eig
#from sympy import Matrix, re, im
#from sympy.physics.quantum import TensorProduct
#from sympy import eye
#from sympy import zeros
#import multiprocessing as mp
import seaborn as sns 
from matplotlib.colors import LinearSegmentedColormap

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

'''def makeList(M, N):
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
    return list'''


    
J = 1.0
h = 0.0
#g = 0.2
#limit = 10e-4

longitudinal_field = sx
transverse_field = sz
N = 7
periodic = False

#dictionary hold site of perturbation as key and perturbing operator as value
#nonHermitianTerms = {0:splus + sminus}

#gamma = []
#realPart = []
#imagPart = []
threshold = np.ones((N,N))
'''for i in range(N):
    for j in range(N):
        threshold[i][j] = np.nan'''
#print(threshold)

for p in range(N):
    for q in range(N):
        if p==q:
            nonHermitianTerms = {p: splus + sminus}
        else:
            nonHermitianTerms = {p: splus, q: sminus}
        thresholdfound = False
        g = 0
        while g < 1.5 and thresholdfound == False :
            ham = findPerturbedHamiltonian(N,J,h,g, periodic, 
                                                     longitudinal_field, transverse_field,
                                                     nonHermitianTerms)
            ev, vecs = np.linalg.eig(ham)
            #print(ev)
            for e in ev:
                if e.imag > 10e-5:
                    threshold[p][q] = g
                    #print(g)
                    thresholdfound = True
            g = g + 0.01    
                
    
    
print(threshold)
'''ax = sns.heatmap(threshold)
ax.invert_yaxis()
plt.show()'''

valuesList = []
for row in threshold:
    for value in row:
        if value not in valuesList:
            valuesList.append(value)
valuesList.sort()            
             
print(valuesList)        


####Plot as grid of dots instead of heatmap
# Get the dimensions of the array
#rows, cols = threshold.shape
'''
# Create a grid of x and y coordinates
x, y = np.meshgrid(range(cols), range(rows))

# Flatten the arrays for plotting
x = x.flatten()
y = y.flatten()
values = threshold.flatten()

# Create the scatter plot
plt.figure(figsize=(6, 6))'''
colors = ["red", "blue", "green", "black"] 

cvals  = [0, 0.25, 0.50,1]

#next 3 lines adapted from stackoverflow
norm=plt.Normalize(min(cvals),max(cvals))
tuples = list(zip(map(norm,cvals), colors))
custom_cmap = LinearSegmentedColormap.from_list("", tuples)
#custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
markers = ["x","o","s","D"]
plt.figure(figsize=(6, 6)) 
for idx in range(len(valuesList)):
    thresholdCopy = threshold.copy()
    
    for row in range(N):
        for column in range(N):
            if thresholdCopy[row][column] > valuesList[idx]*1.01 or thresholdCopy[row][column] < valuesList[idx]*0.99:
                thresholdCopy[row][column] = np.nan   
    rows, cols = thresholdCopy.shape            
    x, y = np.meshgrid(range(cols), range(rows))

    # Flatten the arrays for plotting
    x = x.flatten()
    y = y.flatten()
    values = thresholdCopy.flatten()

    # Create the scatter plot
    #plt.figure(figsize=(6, 6))      
    scatter = plt.scatter(x, y, c=values, s=25, marker = markers[idx], cmap = custom_cmap, norm = norm )      
    
#custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
#scatter = plt.scatter(x, y, c=values, s=25, cmap = custom_cmap)
#scatter.invert_yaxis()
#plt.subplots_adjust(left=-0.1, right=1.1, top=1.1, bottom=-0.1)
plt.margins(0.1)  # 10% margin on all sides

# Add a colorbar to show the scale
#plt.colorbar(scatter, label='Value')

# Invert the y-axis to match the array's row order
#plt.gca().invert_yaxis()

# Add labels and title
#plt.title(r'$\bullet$', color = "tan" + r'$\bullet$')
plt.text(0, 6.25, r'$\diamond$', color = "black", size = 20)
plt.text(1.0, 6.25, r'$\times$', color = "red", size = 10)
plt.text(1.5, 6.25, r'$\gamma_{PT} = 0$')
plt.text(3.0, 6.25, r'$\bullet$', color = "blue", size = 20)
plt.text(3.5, 6.25, r'$\gamma_{PT} = J/4$')
plt.text(5.0, 6.25, r'$\blacksquare$', color = "green", size = 8)
plt.text(5.5, 6.25, r'$\gamma_{PT} = J/2$')
plt.xlabel('p')
plt.ylabel('q')

# Show the plot
plt.show()
#print(gamma_array[0])
end = time.time()
print(end-start)