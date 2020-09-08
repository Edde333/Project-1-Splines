import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.interpolate as inter
import numpy as np
import math


funcx = lambda x: np.cos(x)
funcy = lambda x: np.sin(x)

#get a list of points to fit a spline to as well
N = 10
umin = 0 
umax = 10 
u = np.linspace(umin,umax,N)

#spline fit
xvals = funcx(u)
yvals = funcy(u)

xspline = inter.InterpolatedUnivariateSpline (u, xvals)
yspline = inter.InterpolatedUnivariateSpline (u, yvals)

#figure.subplot.right
mpl.rcParams['figure.subplot.right'] = 0.8

#set up a plot
fig,axes = plt.subplots(1,1,figsize=(9.0,8.0),sharex=True)
ax1 = axes


pind = None #active point
epsilon = 10 #max pixel distance

def update(val):
    global xvals
    global yvals
    global xspline
    global yspline
    # update curve
    l.set_xdata(xvals)
    l.set_ydata(yvals)
    xspline = inter.InterpolatedUnivariateSpline (u, xvals)
    yspline = inter.InterpolatedUnivariateSpline (u, yvals)
    m.set_xdata(xspline(U))
    m.set_ydata(yspline(U))
    # redraw canvas while idle
    print("Update was called")
    fig.canvas.draw_idle()

def reset(event):
    global yvals
    global spline
    #reset the values
    yvals = funcx(u)
    yspline = inter.InterpolatedUnivariateSpline (u, yvals)
    l.set_ydata(yvals)
    m.set_xdata(xspline(U))
    # redraw canvas while idle
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
    xr = np.reshape(xvals,(np.shape(xvals)[0],1))
    yr = np.reshape(yvals,(np.shape(yvals)[0],1))
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
    global yvals
    if pind is None:
        return
    if event.inaxes is None:
        return
    if event.button != 1:
        return
    
    #update yvals
    #print('motion x: {0}; y: {1}'.format(event.xdata,event.ydata))
    xvals[pind] = event.xdata 
    yvals[pind] = event.ydata
    
    # update curve via sliders and draw
    sliders[pind].set_val(yvals[pind])
    fig.canvas.draw_idle()

def onclick(event):
    if event.dblclick:
         print("DOUBLE: ",get_ind_under_point(event) )
         plt.plot([1,2],[2,2])

U = np.linspace(umin,umax,100)
ax1.plot (funcx(U), funcy(U), 'k--', label='original')
l, = ax1.plot (xvals,yvals,color='k',linestyle='none',marker='o',markersize=8)
m, = ax1.plot (xspline(U), yspline(U), 'r-', label='spline')



ax1.set_yscale('linear')
ax1.set_xlim(-3, 3)
ax1.set_ylim(-3,3)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True)
ax1.yaxis.grid(True,which='minor',linestyle='--')
ax1.legend(loc=2,prop={'size':22})

sliders = []

for i in np.arange(N):

    axamp = plt.axes([0.84, 0.8-(i*0.05), 0.12, 0.02])
    # Slider
    s = Slider(axamp, 'p{0}'.format(i), 0, 10, valinit=yvals[i])
    sliders.append(s)

    
for i in np.arange(N):
    #samp.on_changed(update_slider)
    sliders[i].on_changed(update)

axres = plt.axes([0.84, 0.8-((N)*0.05), 0.12, 0.02])
bres = Button(axres, 'Reset')
bres.on_clicked(reset)

fig.canvas.mpl_connect('button_press_event', button_press_callback)
fig.canvas.mpl_connect('button_release_event', button_release_callback)
fig.canvas.mpl_connect('motion_notify_event', motion_notify_callback)
fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

plt.show()

