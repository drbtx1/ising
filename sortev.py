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
metric (inner product is used here) to determine how similar the vectors in each set are to each other'''

def sortEigenvectors(newList, oldList):
    #transpose to simplify indexing
    newListT = columntorow(newList)
    oldListT = columntorow(oldList)
    
    #take advantage of multiplicity - sorted by (degenerate) eigenvalues
    #in other parts of code
    
    N = len(newList)
    counter = 0      
    #calculate multiplicity of degenerate eigenvalues   
    ways = distribution(N)
    for element in ways:
        element = int(element)
        print(element)
        sortedEigenvectorsT = []    
        for eigenvector in range(counter, counter + element):
            bestfit = eigenvector
            print("best = " + str(bestfit))
            bestfitInnerProduct = np.inner(oldListT[eigenvector], newListT[bestfit])
            #could do a swap in newListT instead of appending to another list to avoid double checking
            for eigenvectortocompare in range(counter, counter + element):
                tempinner = np.inner(oldListT[eigenvector], newListT[eigenvectortocompare])
                if tempinner > bestfitInnerProduct:
                    bestfitInnerProduct = tempinner
                    bestfit = eigenvectortocompare
            sortedEigenvectorsT.append(newListT[bestfit]) 
            counter = counter + element
    #transpose to return to column format used by numpy for eigenvalues         
    return(columntorow(sortedEigenvectorsT))        
                
         
            
    
    
    
'''returns a list struct of ways to get n out of N spins up - utilizes symmety
of binomial distribution, prob(N) = prob(0)'''  
def distribution(N):
    W = []
    for n in range(0,N+1):
        W.append(np.math.factorial(N)/(np.math.factorial(N-n)*np.math.factorial(n)))
    return W    

#assumes square list of lists
def columntorow(A):
    length = len(A)
    B = []
    for column in range(0,length):
        temp = []
        for row in range(0,length):
            temp.append(A[row][column])
        B.append(temp)
    print(B)    
    return B    
            
        


newList = [[1,0,0]]
oldList = [[1.01,0,0]]

print(sortEigenvectors(newList,oldList))
    
    
    
    