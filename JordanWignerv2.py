# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 21:37:57 2024

@author: dbtx
"""
import numpy as np

I = np.identity(2)
b = np.array([[0,1],[0,0]])
bdag = np.array([[0,0],[1,0]])

J = 0
h = 1
kappa = 0
L = 3

def find_k_odd(L):
    k = []
    k.append(0)
    multiple = 2*np.pi/L
    n = 1
    while(n*multiple < np.pi):
        k.append(n*multiple)
        k.append(-1*n*multiple)
        n = n + 1
    k.append(n*multiple)
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

#print(find_k_odd(4))

    

def cdag_c(L,first_site, second_site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[first_site] = bdag
    operators[second_site] = b
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    #print(product)    
    return product      
   

def cdag_cdag(L,first_site, second_site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[first_site] = bdag
    operators[second_site] = bdag
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    return product      
  

def c_c(L,first_site, second_site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[first_site] = b
    operators[second_site] = b
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    #negative sign from K 
    return product      
    

def c_cdag(L,first_site, second_site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[first_site] = b
    operators[second_site] = bdag
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    #negative sign from K
    return product      
   
def number(L,site):
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[site] = np.matmul(bdag, b)
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
      
    return product  


        
    

def H_OBC(L):
    sum = 0
    for j in range(0,L-1):
        sum = sum - J*(cdag_c(L,j,j+1) + kappa*cdag_cdag(L,j,j+1) 
                       + c_cdag(L,j,j+1) + kappa*c_c(L,j,j+1)) + h*(2*number(L,j) -1)
    sum = sum + h*(2*number(L,L-1) -1)
    #print(sum)
    return sum

#def H_PBC(L):
#    N = 0
#    for site in range(0,L):
#        N = N + 


print(H_OBC(L))    

values, vecs = np.linalg.eigh(H_OBC(L))
print(vecs)
print(values)
    


        
        