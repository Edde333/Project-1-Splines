#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 11:01:37 2020

@author: oliviaarnfalk
"""
from scipy.sparse import dia_matrix
import numpy as np


def get_diag_matrix(ylen, nbrnodes, stepsize):

    n = nbrnodes
    ex = np.ones(n)
    
    data =  np.array([ex,-4*ex,ex])
    offsets = np.array([-1, 0,  1])
    A = dia_matrix((data, offsets), shape=(n, n)).toarray()
    
    sub1 = np.eye(n, k = -ylen)
    sub2 = np.eye(n, k =ylen)
    
    A = (A + sub1 + sub2)
    
    return A

    

    
    
if __name__ == "__main__":
    
    stepsize = 1
    nbrnodes = 18
    ylen = 3
    A = get_diag_matrix(ylen, nbrnodes, stepsize)