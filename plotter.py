# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 21:10:10 2020

@author: Computer
"""

from matplotlib import pyplot as plt
import numpy as np
from d_spline import d_spline



def get_example_settings():
        
    # Example given
    CONTROL = [(-12.73564, 9.03455),
    (-12.73564, 9.03455),
    (-12.73564, 9.03455),
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
    (39.67575, 17.30712),
    (39.67575, 17.30712),
    (39.67575, 17.30712)]
    
    D = np.array(CONTROL).T
    U_min = 0
    U_max = 1
    U = np.linspace(U_min, U_max, 25*(U_max-U_min)+1)
    U = np.insert(U,0, [0,0])
    U = np.insert(U,-1, [1,1])
    u = np.linspace(U_min, U_max,1000*(U_max-U_min)+1)
    return U,D,u


def update(event):
    # Updates all the graphs and draws the figure again
    global spline_creator
    global U
    global D
    global u
    global curve
    global u_b
    global pind
    global double_clicked

    
    # Update all graphs
    
    # Update D
    D[0,pind] = event.xdata
    D[1,pind] = event.ydata
    
    # Control point graph
    control_graph.set_xdata(D[0])
    control_graph.set_ydata(D[1])
    
    
    # The whole spline
    curve = spline_creator(u,D)
    spline_graph.set_xdata(curve[0])
    spline_graph.set_ydata(curve[1])
    
    # THe intervall points
    intervalls_points =spline_creator(U,D,g=3) 
    intervall_graph.set_xdata(intervalls_points[0])
    intervall_graph.set_ydata(intervalls_points[1])
    
    # Blossom 1
    curve_b1 = spline_creator(u_b, D,g=2) 
    blossom_1.set_xdata(curve_b1[0]) 
    blossom_1.set_ydata(curve_b1[1])
    
    # Blossom 2
    curve_b2 = spline_creator(u_b, D,1) 
    blossom_2.set_xdata(curve_b2[0]) 
    blossom_2.set_ydata(curve_b2[1])
    
    # Blossom 3
    blossom_3.set_xdata(np.array(D[0,double_clicked-2]))
    blossom_3.set_ydata(np.array(D[1,double_clicked-2]))
    
    # Redraw the plots
    fig.canvas.draw_idle()
        
def button_press_callback(event):
    # whenever a mouse button is pressed
    global pind
    if event.inaxes is None:
        return
    if event.button != 1:
        return
    
    # Get whihc of the points (if any) was pressed. (index)
    pind = get_ind_under_point(event) 
    
    # Check if it was a double click
    if event.dblclick:
        global u_b
        global double_clicked
        double_clicked = pind
        # Set new intervall for the blossoms
        u_b = np.linspace(U[double_clicked], U[double_clicked+1], 11)
        update(event)
    

def button_release_callback(event):
    # Whenever a mouse button is released
    global pind
    if event.button != 1:
        return
    pind = None
    
def get_ind_under_point(event):
    # Get the Index of the point pressed (Within the sccuracy given by epsilon)

    # Get distance to points
    tinv = ax1.transData 
    xr = np.reshape(D[0],(np.shape(D[0])[0],1))
    yr = np.reshape(D[1],(np.shape(D[1])[0],1))
    xy_vals = np.append(xr,yr,1)
    xyt = tinv.transform(xy_vals)
    xt, yt = xyt[:, 0], xyt[:, 1]
    d = np.hypot(xt - event.x, yt - event.y)
    indseq, = np.nonzero(d == d.min())
    ind = indseq[0]

    
    if d[ind] >= epsilon:
        ind = None
    
    return ind

def motion_notify_callback(event):
    # Update when things move
    
    if pind is None:
        return
    if event.inaxes is None:
        return
    if event.button != 1:
        return
    
    update(event)
    fig.canvas.draw_idle()
    return 

# Initiate spline variables
U,D,u = get_example_settings()
k = 3 # Order of the algorithm

# Initiate variables for plotting
pind = 0 # selected point
double_clicked = 0
epsilon = 10 #max pixel distance

# Initiate d_spline object
spline_creator = d_spline(k, U)

# Initiate curves

# Whole spline
curve = spline_creator(u,D)
# The interval points
intervalls_points = spline_creator(U,D)
# The first blossom
u_b = np.linspace(U[pind], U[pind+1], 11)
curve_b1 = spline_creator(u_b,D, g = 2) 
# The seconds blossom
curve_b2 = spline_creator(u_b,D, g = 1) 

# print the test example given
print(u[200], ": ", end = "")
print(curve[:,200])


#set up a plot
fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
ax1 = axes


# Initiate Graphs
control_graph, = ax1.plot (D[0],D[1],color='k',linestyle='dotted',marker='o',markersize=8)
spline_graph, = ax1.plot (curve[0], curve[1], 'r-', label='spline')
intervall_graph, = ax1.plot (intervalls_points[0],intervalls_points[1], color='r',linestyle='none',marker='*',markersize=4)
blossom_1,  = ax1.plot (curve_b1[0], curve_b1[1], 'b-', label='Blossom 1')
blossom_2,  = ax1.plot (curve_b2[0], curve_b2[1], 'g-', label='Blossom 2')
blossom_3, =  ax1.plot (np.array(D[0,double_clicked-2]),np.array(D[1,double_clicked -2]), color='c',linestyle='none',marker='s',markersize=10)

# Graph settings
ax1.set_yscale('linear')
ax1.set_xlim(min(D[0])-1, max(D[0])+1)
ax1.set_ylim(min(D[1])-1,max(D[1])+1)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.yaxis.grid(True,which='minor',linestyle='--')
ax1.legend(loc=2,prop={'size':22})

# Configure event callings
fig.canvas.mpl_connect('button_press_event', button_press_callback)
fig.canvas.mpl_connect('button_release_event', button_release_callback)
fig.canvas.mpl_connect('motion_notify_event', motion_notify_callback)
plt.show()

plt.show()














