#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 11:06:46 2020

@author: joel
"""
import numpy as np

class Pollynomial:
    
    
    def __init__(self, vec = np.array([0]) ):
        self.vec = vec
    
    
    def __call__(self, x):
        X = np.array([x**k for k in range(self.order()+1)])
        return self.vec@X
        
    def order(self):
        return self.vec.size - 1
    
    def __str__(self):
        s = "{}".format(self.vec[0])
        for i in range(1,self.order()+1):
            s+= "+{}x^{}".format(self.vec[i], i)
        return s
    
    def mult(self, p2):
        N = self.order()
        M = p2.order()
        mat = np.zeros((M+N+1, M+1))
        for i in range(M+1):
            mat[i:i+N+1,i] = self.vec
        return Pollynomial(mat@p2.vec)
    
    def add(self, p2):
        return Pollynomial(self.vec + p2.vec)
    
    
if __name__ == '__main__':
    v = np.array([1,2,1])
    p1 = Pollynomial(v)
    u = np.array([1,2,1])
    p2 = Pollynomial(u)
    
    print(p1.mult(p2))
    print(p1.add(p2))