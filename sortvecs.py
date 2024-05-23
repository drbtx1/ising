# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:27:09 2024

@author: dabuch
"""


import numpy as np
def sortEigenvectorsasarray(newVectorList, oldVectorList):
    #assumes square arrays. find size of array
    N = len(newVectorList)
    sortedList = np.empty([N,N])
    
    
    #assumes spacing of h is close enough that only nearest neighbor eigenvalues
    #can cross in a single step
    
    
    
    for k in range(0, N):  
            bestfit = 0                      
            bestfitinner = np.inner(squareVector(newVectorList[0:N, 0]),squareVector(oldVectorList[0:N,k]) )
            print(bestfitinner)
            for l in range(1,N):
                tempinner = np.inner((newVectorList[0:N,l]), (oldVectorList[0:N,k]))
                print(tempinner)
                if (tempinner > bestfitinner):
                    bestfit = l
                    bestfitinner = tempinner        
            
                               
                
             
                
            for j in range(N):
                sortedList[j,k] = newVectorList[j,bestfit]   
            #print(sortedList)
                
                
            
    return sortedList    
    #return(columntorow(sortedList))   

def squareVector(originalVector):
    N = len(originalVector)
    squaredVector = np.empty(N)
    for j in range(0,N):
        squaredVector[j] = (originalVector[j])**2
    return squaredVector    

newList = np.array([[1,0,0,0],[0,0,0,0.99],[0,0,1,0],[0,0.98,0,0]])
oldList = np.array([[1.01,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

'''oldList1 = np.array([[-2.28607344e-01,  7.79742596e-02, -7.08677526e-17,  4.93192822e-01,
   1.28248782e-01, -6.19204717e-17,  6.56727441e-01,  5.00680392e-01]
 [ 3.13893156e-01,  1.76959152e-01,  5.00000000e-01, -3.46636041e-01,
  -3.46636041e-01,  5.00000000e-01,  1.76959152e-01,  3.13893156e-01]
 [ 5.00680392e-01, -6.56727441e-01, -1.06656850e-16, -1.28248782e-01,
   4.93192822e-01, -1.48668578e-16,  7.79742596e-02,  2.28607344e-01]
 [-3.13893156e-01,  1.76959152e-01, -5.00000000e-01, -3.46636041e-01,
   3.46636041e-01,  5.00000000e-01, -1.76959152e-01,  3.13893156e-01]
 [ 3.13893156e-01,  1.76959152e-01, -5.00000000e-01, -3.46636041e-01,
  -3.46636041e-01, -5.00000000e-01,  1.76959152e-01,  3.13893156e-01]
 [-5.00680392e-01, -6.56727441e-01,  1.02502178e-17, -1.28248782e-01,
  -4.93192822e-01,  3.58392283e-17, -7.79742596e-02,  2.28607344e-01]
 [-3.13893156e-01,  1.76959152e-01,  5.00000000e-01, -3.46636041e-01,
   3.46636041e-01, -5.00000000e-01, -1.76959152e-01,  3.13893156e-01]
 [ 2.28607344e-01,  7.79742596e-02,  4.86184102e-17,  4.93192822e-01,
  -1.28248782e-01,  1.03515991e-16, -6.56727441e-01,  5.00680392e-01]])

newList2 = np.array([[ 2.30287045e-01, -7.83211203e-02, -3.52979568e-17, -4.95368646e-01,
  -1.27981929e-01, -4.38459260e-17, -6.56192428e-01, -4.98473377e-01]
 [-3.15033466e-01, -1.77872566e-01,  5.00000000e-01,  3.45130795e-01,
   3.45130795e-01,  5.00000000e-01, -1.77872566e-01, -3.15033466e-01]
 [-4.98473377e-01,  6.56192428e-01,  4.08827225e-17,  1.27981929e-01,
  -4.95368646e-01,  5.26841116e-17, -7.83211203e-02, -2.30287045e-01]
 [ 3.15033466e-01, -1.77872566e-01, -5.00000000e-01,  3.45130795e-01,
  -3.45130795e-01,  5.00000000e-01,  1.77872566e-01, -3.15033466e-01]
 [-3.15033466e-01, -1.77872566e-01, -5.00000000e-01,  3.45130795e-01,
   3.45130795e-01, -5.00000000e-01, -1.77872566e-01, -3.15033466e-01]
 [ 4.98473377e-01,  6.56192428e-01, -4.83617687e-17,  1.27981929e-01,
   4.95368646e-01, -1.95390516e-16,  7.83211203e-02, -2.30287045e-01]
 [ 3.15033466e-01, -1.77872566e-01,  5.00000000e-01,  3.45130795e-01,
  -3.45130795e-01, -5.00000000e-01,  1.77872566e-01, -3.15033466e-01]
 [-2.30287045e-01, -7.83211203e-02, -9.86977034e-17, -4.95368646e-01,
   1.27981929e-01, -3.97468634e-17,  6.56192428e-01, -4.98473377e-01]])'''
 
print(sortEigenvectorsasarray(newList,oldList))