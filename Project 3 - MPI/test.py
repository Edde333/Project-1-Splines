# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:19:31 2020

@author: kadde
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

a = np.array([np.linspace(0, 9, 10)])
aext = a
for i in range(9):
    aext = np.concatenate((aext, a))

room2x = np.linspace(0, 1, 10)
room2y = np.linspace(0, 1, 10)

z2 = interpolate.interp2d(room2x, room2y, aext)

znew2 = z2(np.arange(0, 1, 1e-2), np.arange(0, 1, 1e-2))
plt.figure(1)
plt.imshow(znew2)

blackbox = np.zeros((np.shape(znew2)[1], np.shape(znew2)[1]))

a = abs(a-10)
aext = a
for i in range(19):
    aext = np.concatenate((aext, a))

room1x = np.linspace(0, 1, 10)
room1y = np.linspace(0, 2, 20)

z1 = interpolate.interp2d(room1x, room1y, aext)
znew1 = z1(np.arange(0, 1, 1e-2), np.arange(0, 2, 1e-2))
plt.figure(2)
plt.imshow(znew1)    


a = np.array([np.linspace(0,9,10)])
aext = a
for i in range(9):
    aext = np.concatenate((aext,a))

room3x = np.linspace(0,1,np.shape(aext)[1])
room3y = np.linspace(0,1,np.shape(aext)[0])

z3 = interpolate.interp2d(room3x, room3y, aext)
znew3 = z3(np.arange(0,1,1e-2),np.arange(0,1,1e-2))
plt.figure(3)
plt.imshow(znew3)

##### CONCATENATE #####
plot = np.concatenate((blackbox,znew2))
plot = np.concatenate((plot,znew1),1)
plt.figure(4)

right = np.concatenate((znew3,blackbox))
plot = np.concatenate((plot,right),1)
plt.imshow(plot)
