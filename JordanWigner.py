# -*- coding: utf-8 -*-
"""
Created on Wed May 15 18:31:53 2024

@author: dbtx
"""
import numpy as np
import random


length = 7
chain =np.empty(length)
#chain = [1,0,0,1,0,1,1]
for site in range(length):
    chain[site] = random.randint(0,1)
    
#Show original chain for comparison    
print("Original chain: "+ str(chain))   

#variable product is set to 0 if any terms in string of operators disappear
#otherwise gives net phase of expression. I think there would have to be
#product1, product2, etc to use this program to show the result of sums of 
#strings of operators 
product = 1


#If product already equals 0, function b leaves it unchanged. If fermion state
#function acts on is unpopulated, expression becomes 0 (disappears), so "product"
#is set to 0. Otherwise functions populates state with a single particle by
#setting array entry that represents it to 1 
def b(site,temp_chain):
    global product
    if product == 0:
        return temp_chain 
    if temp_chain[site] == 0:
        product = 0
        return temp_chain
    if temp_chain[site] == 1:
       temp_chain[site] = 0
       return temp_chain
    
#If product already equals 0, function bdag leaves it unchanged. If fermion state
#function acts on is populated, expression becomes 0 (disappears), so "product"
#is set to 0. Otherwise functions depopulates state (by removing single particle) by
#setting array entry that represents it to 0         
def bdag(site,temp_chain):
    global product
    if product == 0:
        return temp_chain 
    if temp_chain[site] == 1:
        product = 0
        return temp_chain
    
    if temp_chain[site] == 0:
        temp_chain[site] = 1
        return temp_chain
    

#returns array it took as parameter ("temp_chain) unchanged (for chaining functions.) 
#leaves product equal to 0 if it already has that value. Otherwise multiplies 
#product by net phase of particles to left of site (1 for even number of 
#particles to left, -1 for odd number of particle)
def K(site,temp_chain):
    global product
    if product == 0:
        return temp_chain
    
    for loc in range(0,site):
        if temp_chain[loc] == 1:
            product = -1*product
    return temp_chain

#calls function b, then wraps that in function K. K acts to left of site, so 
#changing value of fermion state at site does not affect result of K, unless b
#causes entire expression to equal 0 
def c(site,temp_chain):
    return (K(site,b(site,chain)))
             
def cdag(site,temp_chain):
    return  (K(site,bdag(site,chain)))             

print("Result:         " + str(cdag(4,c(4,chain))) + "*" + str(product))
       
        
        
        