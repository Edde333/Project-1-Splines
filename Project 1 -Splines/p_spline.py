#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 11:47:16 2020

@author: joel
"""
from pollynomial import Pollynomial
import numpy as np
import matplotlib.pyplot as plt


class p_spline:
    
    
    def __init__(self,U, D, K=3):
        """
        Creates an object which is used to calculate the spline values.
        Parameters
        ----------
        k : Int
            The dimension of the spline (for example 1D, 2D, 3D...)
        u : Array (1xK-4)
            The chosen u-values of the spline
        d : Array (1xK-4)
            The chosen control points of the spline

        Returns
        -------
        None.

        """
        self.K = K
        self.U = U
        self.D = D
        self.dim = 1
        self.P = np.empty((K+1 ,self.dim , U.size-1), dtype = object)

        U_padd = self.expandArray(self.U)
        D_padd = self.expandArray(self.D)
        P_padd = np.empty((K+1 ,self.dim , U.size+2*(K-1) ), dtype = object)

        for i in range(U_padd.size):
            P_padd[0,0,i] = Pollynomial( D_padd[i])
        
        for k in range(1,self.K):
            for i in range(self.K-1-k, self.U.size+)

        

    def expandPollynomial(self,p, k):
        p = np.insert(p, 0, np.array([p[0]  for i in range(k-1)] ) )
        p = np.append(p, np.array(   [p[-1] for i in range(k-1)] ) )
        return p
    
    def expandArray(self, v):
        """
        Expands the vector u by prepending k-1 elements
        equal to the first element in u, and appends 
        k-1 elements equal to the last element in u.
        
        Parameters
        ----------
        u : Array of floats (1,K-4)
            The vector that will be expanded

        Returns
        -------
        v : Array of floats (1,K)
            The extended Array v.

        """
        v = np.insert(v, 0, np.ones(self.K-1)*v[0])
        v = np.append(v, np.ones(self.K-1)*v[-1])
        return v
    

U = np.array([1,2,3])

D = U*2

spline = p_spline(U,D,K= 3)


plt.plot(U,D)
plt.show()



