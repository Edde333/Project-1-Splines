# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
"""
Created on Fri Sep  4 10:47:42 2020

@author: kadde
"""

class d_spline:
    def __init__(self, k, u, d):
        self.k = 3
        self.u = u
        self.d = d
       
    
    def dSpline(self,x):
        self.u = self.expandArray(self.u)
        ret = np.zeros((self.d.ndim,x.shape[0]))
        I = self.findHotInterval(self.u,x)
        # Run dSpline for each row(dimension) of the d matrix.
        if self.d.ndim > 1:
            for dimension in range(self.d.ndim):
                currd = self.expandArray(self.d[dimension])
                left = I
                right = I+1
                depth = 0
                ret[dimension] = self.dRecursive(x, I, currd, depth, left, right)
        else: 
            self.d = self.expandArray(self.d)
            left = I
            right = I+1
            depth = 0
            return self.dRecursive(x, I, self.d, depth, left, right)
        return ret

    def expandArray(self,u):
        u = np.insert(u, 0, np.ones(self.k-1)*u[0])
        u = np.append(u,np.ones(self.k-1)*u[-1])
        return u
    
    def findHotInterval(self,u,x):
        """
        I = np.argmax(u > x) - 1
        if I < k-1: return k-1
        """
        I = np.zeros(len(x), dtype = 'int8')
        for i in range(len(x)):
            index = np.argmax(u >= x[i]) - 1
            if index < self.k-1: index = self.k-1
            I[i] = index
        return I
    
        
        
    def dRecursive(self, x, I, d, depth, left, right):
        if np.any(right-left == self.k+1):
            return d[left+1]
        a = self.alpha(x,left, right)
        d0 = self.dRecursive(x, I, d,  depth + 1, left - 1, right)
        d1 = self.dRecursive(x, I, d, depth + 1, left, right + 1)
        return a*d0 + (1-a)*d1
        
        
        
    def alpha(self, x,left, right):
        alphdivider = (self.u[right]-self.u[left])
        zeroValues = np.where(alphdivider == 0)
        alphdivider[zeroValues] = 1
        alph = (self.u[right] - x)/alphdivider
        alph[zeroValues] = 0
        return alph
    
    
    


# k = 3
# u = np.array([0.5, 2, 5, 8, 9, 11])
# d = np.array([4, 0.5, 3.4, 5.7, 5,2])

# x = np.array([7])
# u = expandArray(u,k)
# d = expandArray(d,k)
# I = findHotInterval(u,x)
# left = I
# right = I+1
# depth = 0

# y = dRecursive(x, I, depth, left, right)
# print(y)

# f = lambda x: dSpline(x)

# x = sp.linspace(1,5,100)

# plt.plot(x,f(x))
# plt.scatter(u,d)
# plt.show()





