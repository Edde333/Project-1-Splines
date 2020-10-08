# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:19:31 2020

@author: kadde
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

a = np.array([np.linspace(0,9,10)])
aext = a
for i in range(10):
    aext = np.concatenate((aext,a))
    
roomx = np.linspace(0,1,10)
roomy = np.linspace(0,1,10)

z = interpolate.interp2d(roomx,roomy,a)
znew = z(np.arange(0,9,1e-2),np.arange(0,9,1e-2))
plt.imshow(znew)    