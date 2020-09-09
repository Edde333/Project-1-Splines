#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:03:40 2020

@author: joel
"""

import numpy as np
import base_functions as bf
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded, inv
from d_spline import d_spline
from math import sin


def get_control_from_inter(U, I):
    L = U.size-3
    eps = np.empty((L+1))
    N = np.empty((L+1), dtype = object)


    for i in range(0,L+1):
        eps[i] = (U[i]+U[i+1]+U[i+2])/3
        N[i]  = bf.getBaseFunc(U, i)
    
    M = np.empty((L+1, L+1))
    
    for i in range(L+1):
        for j in range(L+1):
            M[j,i] = N[i](eps[j])

    
    dx = inv(M)@I[0]
    dx= np.append(dx, [dx[-1], dx[-1]])  
    
    dy = inv(M)@I[1]
    dy= np.append(dy, [dy[-1], dy[-1]])  
    
    D = np.array([dx,dy])
    
    return D



if __name__ == '__main__':
    U = np.linspace(0,10,21)
    U = np.insert(U, 0, np.array([U[0],U[0]]) )
    U = np.insert(U, -1, np.array([U[-1],U[-1]]) )
    U_min = min(U)
    U_max = max(U)
    u = np.linspace(U_min, U_max, 100*int(U_max-U_min)+1)
    
    ix= U*1
    ix = ix[:-2]
    
    iy = np.sin(U)
    iy = iy[:-2]
    
    I = np.array([ix,iy])
    
    D = get_control_from_inter(U, I)
    
    
    ds = d_spline(3,U)
    
    curve = ds(u,D)
    
    plt.scatter(ix,iy)
    plt.plot(curve[0], curve[1])
        
    
    plt.show()
    















