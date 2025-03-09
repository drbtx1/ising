# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 19:57:12 2025

@author: dbtx
"""

import numpy as np
import scipy as sp

c = np.array([[0,1],[0,0]])
cplus = np.array([[0,0],[1,0]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

periodic = False

n = 5
Jplus = np.ones(2**n)
Jminus = np.ones(2**n)
JLplus = 1
JLminus = 1
J = 1
h = 1

def fermionWeights(n):
    weights = np.zeros(n, dtype = int)
    for site in range(0,n):
        weights[site] = 2**(n -1 - site)
    return weights

def ccMatrix(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = c
    operators[j+1] = c
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return product 
def ccMatrixPBC(n):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[0] = c
    operators[n-1] = c
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return -1*product  #account for anticommution   
    


def cpluscMatrix(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = cplus
    operators[j+1] = c
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
        #print(product)
    return product 
def cpluscMatrixPBC(n):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[0] = c
    operators[n-1] = cplus
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return -1*product #account for anti-commution 
 

def cpluscplusMatrix(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = cplus
    operators[j+1] = cplus
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return product  
def cpluscplusMatrixPBC(n):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[0] = cplus
    operators[n-1] = cplus
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return -1*product  

def ccplusMatrix(n, j):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[j] = c
    operators[j+1] = cplus
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return product 

def ccplusMatrixPBC(n):
    operators = []
    for index in range(n):
        operators.append((np.eye(2)))
    operators[0] = cplus
    operators[n-1] = c
    product = operators[0]
    for index in range(1,n):
        product = np.kron(product, operators[index])
    return -1*product     

def hTerm(n):
    operator = np.eye(2**n)
    count = fermionWeights(n)
    #print(count)
    for index in range(2**n):
        term = index
        sum = 0
        for item in count:
            sum += int(term/item)
            term = term % item
        operator[index][index] = 2*sum - n            
    return operator

def HOBC25(n):
    H = np.zeros(2**n)
    for index in range(0,n-1):
        H = H + J*(cpluscMatrix(n,index) 
                              + ccplusMatrix(n,index)) + J*(cpluscplusMatrix(n,index) + ccMatrix(n,index))
    H = -1*H + h*hTerm(n)
    return H

def parityClass(k):
    weights = fermionWeights(n)
    sum = 0
    term = k
    for index in range(len(weights)):
        if int(term/(weights[index])) != 0:
            sum += 1
        term = term % weights[index]        
    return sum % 2     

def N(n):
    sum = np.zeros((2**n,2**n))
    for j in range(n):
        operators = []
        for index in range(n):
            operators.append((np.eye(2)))
        operators[j] = np.array([[0,0],[0,1]])
        product = operators[0]
        for index in range(1,n):
            product = np.kron(product, operators[index])
        sum = sum + product 
    return sum    

#sign check on PBC term: + in paper, c1cN = -cNc1, neither operator changes particle 
#count to left for operator that follows, - overall        
def HPBC25(n):
    H = HOBC25(n) -  J*sp.linalg.expm(1j*np.pi*N(n))@(cpluscMatrixPBC(n) + ccplusMatrixPBC(n)) 
    - J*sp.linalg.expm(1j*np.pi*N(n))@(cpluscplusMatrixPBC(n) +ccMatrixPBC(n))
    return H

def P0(n):
    return (np.real(np.eye(2**n) + sp.linalg.expm(1j*np.pi*N(n))))/2

def P1(n):
    return (np.real(np.eye(2**n) - sp.linalg.expm(1j*np.pi*N(n))))/2

values, vecs = np.linalg.eig(HOBC25(5))
#print(values)
#print(HOBC25(4))
def Hodd(n):
    H = P1(n)@HPBC25(n)@P1(n)
    size = len(H[0,:])
    rows_deleted = 0
    for row in range(size):
        if parityClass(row) == 0:
            H = np.delete(H, row - rows_deleted, 0 )
            #now do column
            H = np.delete(H, row - rows_deleted, 1 )
            rows_deleted += 1
    return H     
def Heven(n):
    H = P0(n)@HPBC25(n)@P0(n)
    size = len(H[0,:])
    rows_deleted = 0
    for row in range(size):
        if parityClass(row) == 1:
            H = np.delete(H, row - rows_deleted, 0 )
            #now do column
            H = np.delete(H, row - rows_deleted, 1 )
            rows_deleted += 1
    return H  
#print(HPBC25(2))
#print(Hodd(2))  
#print(Heven(2))

def HfullPBC25(n):
    H = np.zeros((2**n, 2**n)) #.astype('complex')
    even = Heven(n)
    odd = Hodd(n)
    for row in range(0, 2**(n-1)):
        for column in range(0, 2**(n-1)):
            H[row][column] = even[row][column]
    for row in range(2**(n-1), 2**n):
        for column in range(2**(n-1), 2**n):
            H[row][column] = odd[row- 2**(n-1)][column - 2**(n-1)] 
    return H

#print(HfullPBC25(4))


def find_k_odd(L):
    narray = []
    n = -L/2 + 1
    while (n < L/2):
        narray.append(n)
        n += 1/2
    narray.append(n)    
    k = []
    multiple = 2*np.pi/L
    for idx in range(len(narray)):
        k.append(multiple*narray[idx])
    k.sort()
    return k



def find_k_even(L):
    k = []
    #multiple = np.pi/L
    n = 1
    while((2*n - 1) <= (L - 1) ):
        k.append((2*n - 1)*np.pi/L)
        k.append(-1*(2*n - 1)*np.pi/L)
        n = n + 1
    #k.append(n*multiple)
    k.sort()
    return k        
            


phi = 0
kappa = 0.5
def ck(n,k):
    sum = 0
    for j in range(n):
        operators = []
        for index in range(n):
            operators.append((np.eye(2)))
        operators[j] = c
        product = operators[0]
        for index in range(1,n):
            product = np.kron(product, operators[index])
        sum += np.exp(-1j*k*j)*product
    return np.exp(-1j*phi)*sum/np.sqrt(n)

def cplusk(n,k):
    sum = 0
    for j in range(n):
        operators = []
        for index in range(n):
            operators.append((np.eye(2)))
        operators[j] = cplus
        product = operators[0]
        for index in range(1,n):
            product = np.kron(product, operators[index])
        sum += np.exp(1j*k*j)*product
    return np.exp(1j*phi)*sum/np.sqrt(n)

def nk0(n):
    return cplusk(n,0)@ck(n,0)

def nkpi(n):
    return cplusk(n,np.pi)@ck(n,np.pi)



def Kpluseven(n):
    k = []
    idx = 1
    while idx <= n/2 :
        k.append((2*idx-1)*np.pi/n)
        idx += 1
    return k


def Kplusodd(n):
    k = np.array(find_k_odd(n))
    k =k[k > 0]
    k = k[:-1]
    return k
            

#print(Kpluseven(10))
def Hk0pi(n):
    return -2*J*(nk0(n) - nkpi(n)) + 2*h*(nk0(n) + nkpi(n) - np.eye(2**n))
#print(Hk0pi(3))    

def Hkspace(n,k):
    return 2*(h - J*np.cos(k))*(cplusk(n,k)@ck(n,k) - ck(n,-k)@cplusk(n,-k))\
    -2*kappa*J*np.sin(k)*(1j*np.exp(-2*1j*phi)*cplusk(n,k)@cplusk(n,-k)\
    - 1j*np.exp(2*1j*phi)*ck(n, -k)@ck(n,k))
        
def Hkspaceeven(n):
    sum = np.zeros((2**n, 2**n))
    for k in Kpluseven(n):
        sum += Hkspace(n,k)    
    return sum    
    
def Hkspaceodd(n):
    sum = Hk0pi(n)
    for k in Kplusodd(n):
        sum += Hkspace(n,k)   
    return sum

#print(Hkspaceodd(3))    

def H52(k):
    return 2*kappa*J*(np.sin(2*phi)*np.sin(k)*sx + np.cos(2*phi)*np.sin(k)*sy  \
    + (h - J*np.cos(k))*sz)
        
print(H52(2*np.pi/3)) 

vals, vecs = np.linalg.eig(H52(2*np.pi/3))       
print(vals)
print(vecs)
        
        
    
    
    


