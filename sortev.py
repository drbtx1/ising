# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:02:09 2024

@author: dbtx
"""

import numpy as np


'''takes a list of unsorted eigenvectors (as a list of lists, where the columns 
form the eigenvectors) and the same structure as a list of sorted eigenvectors, 
and returns the list sorted in the same order as the presorted list. This 
assumes that two lists are only incrementally different and there is some  
metric to determine how similar the vectors in each set are to each other'''

def sortEigenvectors(newList, oldList):
    #transpose to simplify indexing
    newListT = newList.T
    oldListT = oldList.T
    
    #take advantage of multiplicity - sorted by (degenerate) eigenvalues
    #in other parts of code
    
    N = len(newList)
    counter = 0             
    
    for element in distribution(N):
        subset = newListT[counter : counter + distribution[element]]
        counter = counter + distribution[element]
        for eigenvector in subset:
            
    
    
    
'''returns a list struct of ways to get n out of N spins up - utilizes symmety
of binomial distribution, prob(N) = prob(0)'''  
def distribution(N):
    W = []
    for n in range(0,N+1):
        print(n)
        W.append(np.math.factorial(N)/(np.math.factorial(N-n)*np.math.factorial(n)))
    return W    

print(distribution(4))
    
    
    
    