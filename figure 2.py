# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 23:05:40 2025

@author: dbtx
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d
import time

start = time.time()


sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

#periodic = False


J = 1
phi = 0
kappa = 1
dh = 0.001
h_list = []
for h in np.arange(0,2,dh):
    h_list.append(h)
k_list = []    
for k in np.arange(-np.pi,np.pi,dh*np.pi):
    k_list.append(k)
    
def H52(k,h,J,phi,kappa):
    return 2*kappa*J*(np.sin(2*phi)*np.sin(k)*sx + np.cos(2*phi)*np.sin(k)*sy + (h - J*np.cos(k))*sz)
        
epsilon_plus = np.zeros((len(h_list),len(k_list)))
epsilon_minus = np.zeros((len(h_list),len(k_list)))   
for k_idx in range(len(k_list)):  
    for h_idx in range(len(h_list)):
        vals, vecs = np.linalg.eigh(H52(h_list[h_idx], k_list[k_idx],J,phi,kappa))
        epsilon_plus[h_idx][k_idx] = max(vals)
        epsilon_minus[h_idx][k_idx] = min(vals)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
#ax = plt.axes(projection ='3d')
ax.plot_surface(h_list, k_list, epsilon_plus, cmap ="viridis" )
ax.plot_surface(h_list, k_list, epsilon_minus, cmap = "plasma")
#fig_size = plt.rcParams["figure.figsize"] #Get current size
#print("Current size:", fig_size) 
ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([0.25, 0.5, 1.25, 1]))

'''fig = plt.figure() 
ax = Axes3D(fig) 
ax.plot(k_list, h_list, epsilon_plus)'''
ax.view_init(elev = 10, azim = 200)
ax.set_title('Figure 2')
plt.show()
        
print(time.time() -start)
