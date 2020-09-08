# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 16:39:12 2020

@author: kadde
"""

from d_spline import *
from base_functions import *
import numpy as np 
import matplotlib.pyplot as plt
import unittest


def testL2Norm(k = 3, u = np.linspace(1,10,10), x = linspace(0,10,100),
               d = np.array([[0, 1, 3, 3, 4, 5, 6, 7, 8, 9, 10],
              [2, 4, 2, 7, 4, 8, 5, 3, 4, 6, 0]])):
    k = 3
    u = np.linspace(1,10,10)
    x = linspace(1,10,100)
    d = np.array([[1, 3, 3, 4, 5, 6, 7, 8, 9, 10],
                  [4, 2, 7, 4, 8, 5, 3, 4, 6, 0]])
    
    
    ds = d_spline(k, u, d)
    res = ds.dSpline(x,3)
    
    u2 = np.insert(u, 0, np.ones(k-1)*u[0])
    u2 = np.append(u2,np.ones(k-1)*u[-1])
    d2x = np.insert(d[0], 0, np.ones(k-1)*d[0][0])
    d2y = np.insert(d[1], 0, np.ones(k-1)*d[1][0])
    N_i = []
    for i in range(len(u2)-2):
        N_i.append(getBaseFunc(u2,i))
        
    base_x = N_i[0](x)*d2x[0]
    base_y = N_i[0](x)*d2y[0]
    
    for i in range(len(N_i)-1):
        base_x +=  N_i[i+1](x)*d2x[i+1]
        base_y +=  N_i[i+1](x)*d2y[i+1]
    
    plt.figure()
    plt.plot(res[0],res[1])
    plt.show()
    plt.scatter(d[0],d[1])
    plt.figure()
    plt.plot(base_x,base_y)
    plt.show()
    plt.scatter(d2x,d2y)
    
    L2norm = max(sqrt((res[0]-base_x)**2 + (res[1] - base_y)**2))
        

    
if __name__ == "__main__":
    testL2Norm()