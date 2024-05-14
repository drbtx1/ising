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

def sortEigenvectors2(newList, oldList):
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
            bestfitInnerProduct = np.inner(oldListT[eigenvector], newListT[eigenvector])
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

def sortEigenvectorsbyDegen(newVectorList, oldVectorList):
    #assumes square "arrays" made of list of lists. find size of array
    N = len(newVectorList)
    sortedList = []
    #find degeneracy, assuming eigenvalues follow binomial distribution
    
    
    #ways = distribution(N)
    ways = [1,2]
    #set counter to index of first element of first degenerate group (composed 
    #of first element only in this distribution)
    
    
    counter = 0
    #step througheach group of degenerate eigenvalues, and sort within that 
    #group so each eigenvector corresponds to closest eigenvector in previous 
    #step of external parameter
    for degen in ways:
        degen = int(degen) 
        print(degen)
        for index in range(counter, counter + degen):
            #convert column vectors to lists
            oldvec = []
            firstvec = []
            for j in range(0,N):
                oldvec.append(oldVectorList[j][index])
                firstvec.append(newVectorList[j][counter])
            #print(oldvec)
            #print(firstvec)
            bestfit = counter
            bestfitinner = np.inner(firstvec,oldvec)
            for vectortocompare in range(counter + 1, counter + degen):
                nextvec = []
                for j in range(0,N):
                    nextvec.append(newVectorList[j][vectortocompare])
                tempinner = np.inner(nextvec,oldvec)
                if tempinner > bestfitinner:
                    bestfit = vectortocompare
                    bestfitinner = tempinner
            tempsorted = []
            for j in range(0,N):
                tempsorted.append(newVectorList[j][bestfit])
            sortedList.append(tempsorted)
        counter = counter + degen
    return(np.array(sortedList).T)    


def sortEigenvectorsbypairs(newVectorList, oldVectorList):
    #assumes square "arrays" made of list of lists. find size of array
    N = len(newVectorList)
    sortedList = []
    
    #assumes spacing of h is close enough that only nearest neighbor eigenvalues
    #can cross in a single step
    
    
    
    for index in range(0, N):
            #print(index)
            #convert column vectors to lists
            oldvec = []
            firstvec = []
            for j in range(0,N):
                oldvec.append(oldVectorList[j][index])
                firstvec.append(newVectorList[j][index])
            #print("oldvec" + str(oldvec))
            #print("firstvec" +str(firstvec))
            bestfit = index
            #print("bestfit" + str(bestfit) )
            bestfitinner = np.inner(firstvec,oldvec)
            #print("bestfitinner" + str(bestfitinner) )
            
            if index == 0:
                vectortocompare = 1 
                nextvec = []
                for j in range(0,N):
                    nextvec.append(newVectorList[j][vectortocompare])
                #print("nextvec" + str(nextvec))    
                tempinner = np.inner(nextvec,oldvec)
                if tempinner > bestfitinner:
                    bestfit = vectortocompare
                    bestfitinner = tempinner
                tempsorted = []
                for j in range(0,N):
                    tempsorted.append(newVectorList[j][bestfit])
                sortedList.append(tempsorted)
                #print(sortedList)
                
            elif index == (N - 1):
                vectortocompare = N -2 
                prevvec = []
                for j in range(0,N):
                    prevvec.append(newVectorList[j][vectortocompare])
                tempinner = np.inner(prevvec,oldvec)
                if tempinner > bestfitinner:
                    bestfit = vectortocompare
                    bestfitinner = tempinner
                tempsorted = []
                for j in range(0,N):
                    tempsorted.append(newVectorList[j][bestfit])
                sortedList.append(tempsorted)
                #print(sortedList)    
            else: 
                for vectortocompare in (index - 1, index + 1):
                    nextvec = []
                    for j in range(0,N):
                        nextvec.append(newVectorList[j][vectortocompare])
                    tempinner = np.inner(nextvec,oldvec)
                    if tempinner > bestfitinner:
                        bestfit = vectortocompare
                        bestfitinner = tempinner
                tempsorted = []
                for j in range(0,N):
                    tempsorted.append(newVectorList[j][bestfit])
                sortedList.append(tempsorted)   
                #print(sortedList)
                
                
            
        
    return(np.array(sortedList).T)             
            
def sortEigenvectorscompareall(newVectorList, oldVectorList):
    #assumes square "arrays" made of list of lists. find size of array
    N = len(newVectorList)
    sortedList = []
    
    #assumes spacing of h is close enough that only nearest neighbor eigenvalues
    #can cross in a single step
    
    
    
    for index in range(0, N):
            #print(index)
            #convert column vectors to lists
            oldvec = []
            firstvec = []
            for j in range(0,N):
                oldvec.append(oldVectorList[j][index])
                firstvec.append(newVectorList[j][index])
            #print("oldvec" + str(oldvec))
            #print("firstvec" +str(firstvec))
            bestfit = index
            #print("bestfit" + str(bestfit) )
            bestfitinner = np.inner(firstvec,oldvec)
            #print("bestfitinner" + str(bestfitinner) )
            
            
                
               
             
            for vectortocompare in range(0, N):
                if vectortocompare == index:
                    continue
                #print(vectortocompare)
                nextvec = []
                for j in range(0,N):
                    nextvec.append(newVectorList[j][vectortocompare])
                tempinner = np.inner(nextvec,oldvec)
                if tempinner > bestfitinner:
                    bestfit = vectortocompare
                    bestfitinner = tempinner
                #print(bestfit)    
            '''for vectortocompare in (index + 1, N ):
                nextvec = []
                for j in range(0,N):
                    nextvec.append(newVectorList[j][vectortocompare])
                tempinner = np.inner(nextvec,oldvec)
                if tempinner > bestfitinner:
                    bestfit = vectortocompare
                    bestfitinner = tempinner'''        
            tempsorted = []
            for j in range(0,N):
                tempsorted.append(newVectorList[j][bestfit])
            sortedList.append(tempsorted)   
            #print(sortedList)
                
                
            
    return np.array(sortedList).T    
    #return(columntorow(sortedList))             
            
    
    
    
'''returns a list struct of ways to get n out of N spins up - utilizes symmety
of binomial distribution, prob(N) = prob(0)'''  
def distribution(N):
    W = []
    for n in range(0,N+1):
        W.append(np.math.factorial(N)/(np.math.factorial(N-n)*np.math.factorial(n)))
    return W    

#assumes square list of lists
def columntorow(A):
    length = len(A[0])
    height = len(A)
    
    B = []
    for column in range(0,length ):
        temp = []
        for row in range(0,height):
            temp.append(A[row][column])
        B.append(temp)
    
    return B  

 
         
            
        


newList = [[1,0,0],[0,0,0.99],[0,1,0]]
oldList = [[1.01,0,0],[0,1,0],[0,0,1]]
 
#print(sortEigenvectorscompareall(newList,oldList))
    
#print(columntorow([[1,2,3],[4,5,6],[7,8,9]]))    
    
    