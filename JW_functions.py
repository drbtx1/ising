# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 08:00:33 2025

@author: dbtx
"""
import numpy as np

periodic = False


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
#print(numberOperator([1,1,1], 1))
def cpluscplus(fermionChain, j):
    if j != len(fermionChain) - 1:
        fermionChain[j] = fermionChain[j]^1
        fermionChain[j+1] = fermionChain[j+1]^1
    else:
        fermionChain[j] = fermionChain[j]^1
        fermionChain[0] = fermionChain[0]^1
        fermionChain = -1*fermionChain
    
    return(fermionChain)
#print(cpluscplus([0,0,0],0))    

def cplusc(fermionChain, j):
    fermionChain[j] = fermionChain[j]^1
    fermionChain[j+1] = 0
    return(fermionChain)
    
def cc(fermionChain, j):
    fermionChain[j] = 0
    fermionChain[j+1] = 0
    return(-1*fermionChain)

def ccplus(fermionChain, j):
    fermionChain[j] = 0
    fermionChain[j+1] = fermionChain[j+1]^1
    return(-1*fermionChain)

        
                 
        
        
    