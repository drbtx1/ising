# -*- coding: utf-8 -*-
"""
Created on Mon May 13 10:42:28 2024

@author: dbtx
"""
import numpy as np
def sortEigenvectorsasarray(newVectorList, oldVectorList):
    #assumes square "arrays" made of list of lists. find size of array
    N = len(newVectorList)
    sortedList = np.empty([N,N])
    
    #assumes spacing of h is close enough that only nearest neighbor eigenvalues
    #can cross in a single step
    
    
    
    for index in range(0, N):
            
            bestfit = index
            #print("bestfit" + str(bestfit) )
            bestfitinner = np.inner(newVectorList.T[index],oldVectorList.T[index])
            #print("bestfitinner" + str(bestfitinner) )
                
                
               
             
            for vectortocompare in range(0, N):
                if vectortocompare == index:
                    continue
                #print(vectortocompare)
               
                tempinner = np.inner(newVectorList.T[vectortocompare],oldVectorList.T[index])
                if tempinner > bestfitinner:
                    bestfit = vectortocompare
                    bestfitinner = tempinner
                #print(bestfit)    
                  
            for j in range(N):
                sortedList[j,index] = newVectorList[j,bestfit]   
            #print(sortedList)
                
                
            
    return sortedList    
    #return(columntorow(sortedList))   



newList = np.array([[1,0,0,0],[0,0,0,0.99],[0,0,1,0],[0,0.98,0,0]])
oldList = np.array([[1.01,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
 
print(sortEigenvectorsasarray(newList,oldList))          