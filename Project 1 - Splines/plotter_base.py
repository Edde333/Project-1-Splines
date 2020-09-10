# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 21:10:10 2020

@author: Computer
"""

import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.interpolate as inter
import numpy as np
import math
from d_spline import d_spline
import base_functions as bf



# Example given
CONTROL = [(-12.73564, 9.03455),
(-26.77725, 15.89208),
(-42.12487, 20.57261),
(-15.34799, 4.57169),
(-31.72987, 6.85753),
(-49.14568, 6.85754),
(-38.09753, -1e-05),
(-67.92234, -11.10268),
(-89.47453, -33.30804),
(-21.44344, -22.31416),
(-32.16513, -53.33632),
(-32.16511, -93.06657),
(-2e-05, -39.83887),
(10.72167, -70.86103),
(32.16511, -93.06658),
(21.55219, -22.31397),
(51.377, -33.47106),
(89.47453, -33.47131),
(15.89191, 0.00025),
(30.9676, 1.95954),
(45.22709, 5.87789),
(14.36797, 3.91883),
(27.59321, 9.68786),
(39.67575, 17.30712)]

#D = np.array(CONTROL)
D = CONTROL

U_min = 0
U_max = 1
U = np.linspace(U_min, U_max, 25*(U_max-U_min)+1)



pind = 0 # selected point
epsilon = 10 #max pixel distance

"""
U = np.array([1.,2,3,4,5])
D_x = np.array([1.,2,3,4,5])
D_y = np.array([1.,2,3,4,5])
"""
k = 3 # Order of the algorithm

u = np.linspace(U_min, U_max,1000*(U_max-U_min)+1)

# The actual spline
curve = bf.getCubicSpline(u,U,D)

# print the test example given
print(u[200], ": ", end = "")
print(curve[:,200])

#figure.subplot.right
mpl.rcParams['figure.subplot.right'] = 0.8
#set up a plot
fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
ax1 = axes


def update(val):
    # Updates all the graphs and draws the figure again
    global U
    global D
    global u
    global curve
    
    
    points.set_xdata(D[:,0])
    points.set_ydata(D[:,1])
    
    
    curve = bf.getCubicSpline(u,U,D)
    graph.set_xdata(curve[0])
    graph.set_ydata(curve[1])
    


    spline = d_spline(k, U, D )
    intervalls_points =spline.dSpline(U,3) 
    intervalls.set_xdata(intervalls_points[0])
    intervalls.set_ydata(intervalls_points[1])

    spline = d_spline(k, U, D )
    curve_b1 = spline.dSpline(u_b1,2) 
    blossom_1.set_xdata(curve_b1[0]) 
    blossom_1.set_ydata(curve_b1[1])
    
    spline = d_spline(k, U, D )
    curve_b2 = spline.dSpline(u_b2,1) 
    blossom_2.set_xdata(curve_b2[0]) 
    blossom_2.set_ydata(curve_b2[1])

    fig.canvas.draw_idle()
        
def button_press_callback(event):
    'whenever a mouse button is pressed'
    global pind
    if event.inaxes is None:
        return
    if event.button != 1:
        return
    #print(pind)

    pind = get_ind_under_point(event) 
    if event.dblclick:
        global u_b1
        global u_b2
        u_b1 = np.linspace(U[pind], U[pind+1], 11)
        u_b2 = np.linspace(U[pind], U[pind+1], 11)
        update(event)
    print("Pressed thing")       
    

def button_release_callback(event):
    'whenever a mouse button is released'
    global pind
    if event.button != 1:
        return
    pind = None
    
def get_ind_under_point(event):
    'get the index of the vertex under point if within epsilon tolerance'

    # display coords
    #print('display x is: {0}; display y is: {1}'.format(event.x,event.y))
    t = ax1.transData.inverted()
    tinv = ax1.transData 
    xy = t.transform([event.x,event.y])
    #print('data x is: {0}; data y is: {1}'.format(xy[0],xy[1]))
    xr = np.reshape(D[0],(np.shape(D[0])[0],1))
    yr = np.reshape(D[1],(np.shape(D[1])[0],1))
    xy_vals = np.append(xr,yr,1)
    xyt = tinv.transform(xy_vals)
    xt, yt = xyt[:, 0], xyt[:, 1]
    d = np.hypot(xt - event.x, yt - event.y)
    indseq, = np.nonzero(d == d.min())
    ind = indseq[0]

    #print(d[ind])
    if d[ind] >= epsilon:
        ind = None
    
    #print(ind)
    return ind

def motion_notify_callback(event):
    'on mouse movement'

    global D
    if pind is None:
        return
    if event.inaxes is None:
        return
    if event.button != 1:
        return
    
    #update yvals
    #print('motion x: {0}; y: {1}'.format(event.xdata,event.ydata))
    #D_x[pind] = event.xdata
    #D_y[pind] = event.ydata
    D[0,pind] = event.xdata
    D[1,pind] = event.ydata
    # update curve via sliders and draw
    #sliders[pind].set_val(D_y[pind])
    update(event)
    fig.canvas.draw_idle()


def onclick(event):
    if event.dblclick:
        global u_b1
        #u_b1 = np.linspace(U[pind], U[pind+1], 11)



points, = ax1.plot (D[0],D[1],color='k',linestyle='none',marker='o',markersize=8)
graph, = ax1.plot (curve[0], curve[1], 'r-', label='spline')
intervalls, = ax1.plot (intervalls_points[0],intervalls_points[1], color='r',linestyle='none',marker='*',markersize=4)
blossom_1,  = ax1.plot (curve_b1[0], curve_b1[1], 'b-', label='spline')
blossom_2,  = ax1.plot (curve_b2[0], curve_b2[1], 'g-', label='spline')



ax1.set_yscale('linear')
ax1.set_xlim(min(D[0])-1, max(D[0])+1)
ax1.set_ylim(min(D[1])-1,max(D[1])+1)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.yaxis.grid(True,which='minor',linestyle='--')
ax1.legend(loc=2,prop={'size':22})




fig.canvas.mpl_connect('button_press_event', button_press_callback)
fig.canvas.mpl_connect('button_release_event', button_release_callback)
fig.canvas.mpl_connect('motion_notify_event', motion_notify_callback)
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

plt.show()














