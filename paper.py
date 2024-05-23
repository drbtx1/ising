# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:12:29 2024

@author: dabuch
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
#from sortev import *
#from sortvecs import *


#Following lines make plots look a little more Latex-like
#mpl.rcParams['mathtext.fontstyle'] = 'cm' 
mpl.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmss10'
mpl.rcParams['font.monospace'] = 'cmtt10'
mpl.rcParams["axes.formatter.use_mathtext"] = True # to fix the minus signs


#Some shorthand for frequently used matrices (splus, sminus not used here)
I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
splus = (sx + 1j*sy)/2
sminus = (sx - 1j*sy)/2


#First three  functions are used to build Hamiltonian and find eigenvalues 
#and eigenvectors 

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
  
    

                
#takes number of sites N, energy scaling J and h, and boolean periodic
#(boundary conditions)
#return lists of eigenvalues and eigenfunctions for corresponding Hamiltonian 
def findHamiltonian(N,J,h, periodic,longitudinal_field,transverse_field):   
    tensor_sum_J = build_J(N, periodic,longitudinal_field)
    tensor_sum_h = build_h(N, periodic, transverse_field)
    H = -J*tensor_sum_J - h*tensor_sum_h  
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

#expects a tensor products expectation value operator and eigenvector
#returns an expectation value
def findExpectationValue(pauliOperatorProduct, eigenvector):
    return np.matmul(eigenvector.conj().T,np.matmul(pauliOperatorProduct, eigenvector))     


#calculates mean of expectation across N sites of given eigenvector
def findAndSaveMagnetization(N,periodic,eV,J,maxh,steps,longitudinal_field, transverse_field, direction_of_M):
    h_per_J = []
    operator_list = []
    avgExpValue = []
    

#build list of operators to find expectation value of nth site
#nth element in each list holds tensor product of Pauli operator at nth site
#and all other sites holding identity operator
    for site in range(0,N):
        operator_list.append(expectationOperator(site,N,direction_of_M))
     
    last = []   #used in attempts to sort eigenvectors in same order across iterations
    for step in range(steps): 
        #hold expectation values at nth site 
        expectation_at_site = [] 
        h = maxh*step/steps
        
        eigenvalues, eigenvectors = findHamiltonian(N,J,h,periodic,longitudinal_field, transverse_field)
        #assure sorted in decreasing order
        idx = eigenvalues.argsort()[::-1]   
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:,idx]
        
        #following 2 lines are when eigenvectors are sorted (else commented out)
        #if step != 0:
            #eigenvectors = sortEigenvectorsasarray(eigenvectors,last)
            
        last = eigenvectors    #also used to attempt to sort
        
        for site in range(N):
            expectation_at_site.append(findExpectationValue(operator_list[site], eigenvectors[:,eV]))
            
        avg = sum(expectation_at_site)/len(expectation_at_site)   
        
        h_per_J.append(h/J)
        avgExpValue.append(avg)  #average across N sites, all with same eigenvector
       
        
        
        
    #Following lines were used to save raw data 
    #title = "N=" +str(N)+",eigenvector = " + str(eV) + ", periodic =" + str(periodic)
    #filetype = '.csv'
    #dest = 'C:/Users/dabuch/Ising Data/Magnetization/PeriodicTrue/JxhzMx/N' + str(N) + '/' + title + filetype
    #data = {"h/J": h_per_J, "Expectation Value averaged over all sites": avgExpValue}
    #df = pd.DataFrame(data)
    #df.to_csv(dest)'''
    
    plt.xlabel("h/J")
    plt.ylabel("M")
    #plt.title(title)   
    plt.plot(h_per_J, avgExpValue)
    
    
    #Next 3 lines used to save plots
    #image_filetype = '.png'
    #dest = 'C:/Users/dabuch/Ising Data/Magnetization/PeriodicTrue/JxhzMx/N' + str(N) + '/' + title + image_filetype
    #plt.savefig(dest)
    plt.show()

#Next pair of functions were used to try to sort eigenvectors. They would have 
#to be placed before functions that called them to work (I called from another file)

#def sortEigenvectorsasarray(newVectorList, oldVectorList):
    #assumes square "arrays" made of list of lists. find size of array
    #N = len(newVectorList)
    #sortedList = np.empty([N,N])
        
    #for k in range(0, N):  
    #        bestfit = 0                      
    #        bestfitinner = np.inner(squareVector(newVectorList[0:N, 0]),squareVector(oldVectorList[0:N,k]) )
    #        print(bestfitinner)
    #        for l in range(1,N):
    #            tempinner = np.inner(squareVector(newVectorList[0:N,l]), squareVector(oldVectorList[0:N,k]))
    #            print(tempinner)
    #            if (tempinner > bestfitinner):
    #                bestfit = l
    #                bestfitinner = tempinner        
                                 
      
            #for j in range(N):
                #sortedList[j,k] = newVectorList[j,bestfit]   
            #print(sortedList)
                
                
            
    return sortedList.T    
    #return(columntorow(sortedList))   

def squareVector(originalVector):
    N = len(originalVector)
    squaredVector = np.empty(N)
    for j in range(0,N):
        squaredVector[j] = (originalVector[j])**2
    return squaredVector        
            
                    
        

