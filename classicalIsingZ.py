# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 20:34:33 2025

@author: dbtx
"""

import numpy as np

J = 1
h = 0.5
beta = 1
N = 5
periodic = False



def makePermutations(oldPerm, site):
    temp = oldPerm.copy()
    temp[site] = 1
    if site == N-1:        
        
        permutations.append(oldPerm)
        permutations.append(temp)
        return
    else:
        makePermutations(oldPerm, site + 1)
        makePermutations(temp, site + 1)
#makePermutations(allUp,0)       
#print(permutations)

def findHam(perm):
    sum = 0
    for site in range(len(perm) - 1):
        if perm[site] == perm[site + 1]:
            sum += J
        else:
            sum -= J
    if periodic: 
        if perm[0] == perm[len(perm) - 1]:
            sum += J
        else:
            sum -= J
    for site in range(len(perm)): 
        if perm[site] == 0:
            sum += h
        else:
            sum -= h
    return sum    
        
#print(findHam(allUp))

def Z():
    sum = 0
    

    for perm in permutations:
        
        sum += np.exp(-beta*findHam(perm))
        #print(perm)
        #print(sum)
    return sum
allUp = []
for n in range(N):
    allUp.append(0)
permutations = []
makePermutations(allUp,0)
print(Z())  
print(permutations)  
        
        
    