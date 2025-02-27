# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 08:00:33 2025

@author: dbtx
"""
import numpy as np

c = np.array([[0,1],[0,0]])
cplus = np.array([[0,0],[1,0]])

periodic = False

n = 5
Jplus = np.ones(2**n)
Jminus = np.ones(2**n)
h = 0.5

def fermionWeights(n):
    weights = np.zeros(n, dtype = int)
    for site in range(0,n):
        weights[site] = 2**(n -1 - site)
    return weights
#print(fermionWeights(3))

    
def spinToFermion(spinChain):
    N = len(spinChain)
    n = int(np.log2(N))
    fermionChain = np.zeros(n, dtype = int)
    weights = fermionWeights(n)
    for site in range(0,N):
        spinValue = N - 1 - site
        print(spinValue)
        temp = np.zeros(n, dtype = int)
        for fermloc in range(0,n):
            temp[fermloc] = spinValue/weights[fermloc]
            spinValue = spinValue%weights[fermloc]
        fermionChain = fermionChain + spinChain[site]*temp
    return fermionChain    
        
#print(spinToFermion([1,1,0,0,0,0,0,0]/np.sqrt(2)))       

def numberOperator(fermionChain, j):
    return fermionChain[j]

def particleCount(fermionChain):
    sum = 0
    for j in range(len(fermionChain)):
        sum += numberOperator(fermionChain, j) 
    return sum    
#print(numberOperator([1,1,1], 1))

'''
def cpluscplus(fermionChain, j):
    tempchain = np.copy(fermionChain)
    next = j+1
    factor = 0
    if periodic == True and j == len(fermionChain) - 1:
        next = 0
        factor = (particleCount(fermionChain) % 2) + 1
    tempchain[j] = fermionChain[j]^1
    tempchain[next] = fermionChain[next]^1
    return(((-1)**factor)*tempchain)
    
    
#print(cpluscplus([0,0,0],0))    

def cplusc(fermionChain, j):
    tempchain = np.copy(fermionChain)
    next = j+1
    factor = 0
    if periodic == True and j == len(fermionChain) - 1:
        next = 0
        factor = (particleCount(fermionChain) % 2) + 1
    tempchain[j] = (fermionChain[j])^1
    tempchain[next] = 0
    return(((-1)**factor)*tempchain)
    
def cc(fermionChain, j):
    tempchain = np.copy(fermionChain)
    next = j+1
    factor = 0
    if periodic == True and j == len(fermionChain) - 1:
        next = 0
        factor = (particleCount(fermionChain) % 2) + 1
    tempchain[j] = 0
    tempchain[next] = 0
    return(((-1)**factor)*tempchain)

def ccplus(fermionChain, j):
    tempchain = np.copy(fermionChain)
    next = j+1
    factor = 0
    if periodic == True and j == len(fermionChain) - 1:
        next = 0
        factor = (particleCount(fermionChain) % 2) + 1
    tempchain[j] = 0
    tempchain[next] = fermionChain[next]^1
    return(int(((-1)**factor))*tempchain)
#print(particleCount([1,1,1,1])) '''



def ccMatrixOBC(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = c
    operators[j+1] = c
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return product    
      
    


def cpluscMatrixOBC(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = cplus
    operators[j+1] = c
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
        print(product)
    return product  

def cpluscplusMatrixOBC(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = cplus
    operators[j+1] = cplus
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return product  

def ccplusMatrixOBC(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = c
    operators[j+1] = cplus
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return product  

def hTerm(n,j):
    operator = np.eye(2**n)
    operator[2**(j+1) - 1][2**(j+1) - 1] = -1
    return operator

#print(hTerm(4,3))

def HOBC25(n):
    H = np.zeros(2**n)
    for index in range(0,n-1):
        H = H + Jplus[index]*(cpluscMatrixOBC(n,index) - ccplusMatrixOBC(n,index)) 
        + Jminus[index]*(cpluscplusMatrixOBC(n,index) - ccMatrixOBC(n,index))
        + h*hTerm(n, index)
    H = -1*H + h*hTerm(n,n-1)
    return H

"""def HPBC25(n):
    index = n-1
    H = HOBC25(n) - (Jplus[index]*(cpluscMatrixOBC(n,index) - ccplusMatrixOBC(n,index)) 
    + Jminus[index]*(cpluscplusMatrixOBC(n,index) - ccMatrixOBC(n,index)))
    
def HABC25(n):
    index = n-1
    H = HOBC25(n) + (Jplus[index]*(cpluscMatrixOBC(n,index) - ccplusMatrixOBC(n,index)) 
    + Jminus[index]*(cpluscplusMatrixOBC(n,index) - ccMatrixOBC(n,index)))  """  
    
    
def parityClass(k):
    weights = fermionWeights(n)
    sum = 0
    term = k
    for index in range(len(weights)):
        if int(term/(weights[index])) != 0:
            sum += 1
        term = term % weights[index]
        
    return sum % 2     
        
def HPBC25(n):
    Heven = np.zeros(2**(n-1))
    Hodd = np.zeros(2**(n-1))
    
    
print(parityClass(6))    
    
    
       
            
            
        
    


        
               
        
        
    