# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 21:37:57 2024

@author: dbtx
"""
import numpy as np

I = np.identity(2)
b = np.array([[0,1],[0,0]])
bdag = np.array([[0,0],[1,0]])

J = 1
h = 0.1
kappa = 0.5
L = 3
phi = 0

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


        
    
###DOUBLE CHECK SIGN ON HC
def H_OBC(L):
    sum = 0
    for j in range(0,L-1):
        sum = sum - J*(cdag_c(L,j,j+1) + kappa*cdag_cdag(L,j,j+1) 
                       - c_cdag(L,j,j+1) - kappa*c_c(L,j,j+1)) + h*(2*number(L,j) -1)
    sum = sum + h*(2*number(L,L-1) -1)
    #print(sum)
    return sum

def H_PBC(L, parity = 1):
    N = parity % 2
    return H_OBC(L) + ((-1)**N)*(- J*(cdag_c(L,L-1,1) + kappa*cdag_cdag(L,L-1,1) 
                   - c_cdag(L,L-1,1) - kappa*c_c(L,L-1,1)))


'''print(H_PBC(L))    

values, vecs = np.linalg.eigh(H_OBC(L))
print(vecs)
print(values)'''

def cdag_cdag_same_site(L,site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[site] = np.matmul(bdag,bdag)
    #operators[second_site] = bdag
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    return product 

def c_c_same_site(L,site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[site] = np.matmul(b,b)
    #operators[second_site] = bdag
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    return product 
#print(c_c_same_site(4,1))

#could just number for following
def cdag_c_same_site(L,site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[site] = np.matmul(bdag,b)
    #operators[second_site] = bdag
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    return product 

def c_cdag_same_site(L,site):    
    operators = []
    for site in range(0,L):
        operators.append(I)
    operators[site] = np.matmul(b,bdag)
    #operators[second_site] = bdag
    product = operators[0]
    for site in range(1,L):
        product = np.kron(product, operators[site])
    return product 

def H_PBC_k_space_odd(L):
    #if parity == 1:
    ks = find_k_odd(L)
    m = ks.index(0)
    n = index(np.pi)
    sum = -2*J*(number(L,m) - number(L,n)) + 2*h*(number(L,m) + number(L,n) + np.identity(2**L) )
    kpos = []
    for k in ks:
        if k>0 and k != np.pi:
            kpos.append(k)
    
    for k in kpos:
        m = ks.index(k)
        n = ks.index(-k)
        sum = 2*(h - J*np.cos(k))*(cdag_c_same_site(L,m) - c_cdag_same_site(L,n)) 
        + 2*kappa*J*np.sin(k)*(1j*(np.cos(-2*phi) + 1j*np.sin(-2*phi))*cdag_cdag(L,m,n) 
       - (1j*(np.cos(2*phi) + np.sin(2*phi)))*c_c(L,m,n))
    return sum    

def H_PBC_k_space_even(L):
    #if parity == 1:
    ks = find_k_even(L)
    
    sum = 0
    kpos = []
    for k in ks:
        if k > 0: 
            kpos.append(k)
    
    for k in kpos:
        m = ks.index(k)
        n = ks.index(-k)
        sum = 2*(h - J*np.cos(k))*(cdag_c_same_site(L,m) - c_cdag_same_site(L,n)) 
        + 2*kappa*J*np.sin(k)*(1j*(np.cos(-2*phi) + 1j*np.sin(-2*phi))*cdag_cdag(L,m,n)) 
        - (1j*(np.cos(2*phi) + np.sin(2*phi))*c_c(L,m,n))        
    return sum
    


print(find_k_odd(4))
print(find_k_even(4))
print(H_PBC_k_space_even(4))      

    
        
        
    


        
        