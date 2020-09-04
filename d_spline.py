# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
"""
Created on Fri Sep  4 10:47:42 2020

@author: kadde
"""
def expandArray(u, k):
    u = np.insert(u, 0, np.ones(k-1)*u[0])
    u = np.append(u,np.ones(k-1)*u[-1])
    return u

def findHotInterval(u,x):
    """
    I = np.argmax(u > x) - 1
    if I < k-1: return k-1
    """
    I = np.zeros(len(x), dtype = 'int8')
    for i in range(len(x)):
        index = np.argmax(u >= x[i]) - 1
        if index < k-1: index = k-1
        I[i] = index
    return I

    
    
def dRecursive(x, I, depth, left, right):
    if np.any(right-left == k+1):
        return d[left+1]
    a = alpha(x,left, right)
    d0 = dRecursive(x, I, depth + 1, left - 1, right)
    d1 = dRecursive(x, I, depth + 1, left, right + 1)
    return a*d0 + (1-a)*d1
    
    
    
def alpha(x,left, right):
    return (u[right] - x)/(u[right]-u[left])


def dSpline(x):
    I = findHotInterval(u,x)
    left = I
    right = I+1
    depth = 0
    y = dRecursive(x, I, depth, left, right)
    return y


k = 3
u = np.array([1, 2, 3, 4, 5, 6, 7])
d = np.array([1, 3, 3,1, 0,0,1])



x = np.array([5])
u = expandArray(u,k)
d = expandArray(d,k)
I = findHotInterval(u,x)
left = I
right = I+1
depth = 0

y = dRecursive(x, I, depth, left, right)
print(y)







f = lambda x: dSpline(x)

x = sp.linspace(1,7,100)

plt.plot(x,f(x))
plt.scatter(u,d)
plt.show()





