import matplotlib.pyplot as plt
import sys
import numpy as np

def plot_newton_2d(function, points, levels = None, axes = None):
    """
    
    Parameters
    ----------
    function : Python function from R^2 to R that takes an array of length
        2 as its argument.
        The function on which you want to find a minimum.
    points : (Nx2)-array
        The recursive points that you wish to plot.
    levels : array-like, optional
        An array with the (increasing) levels that the contour plot
        is drawing. If not specified 50 equally spaced (value-wise)
        contours are drawn. Default is None.
    axes : plot-axes, optional
        The axes on which you want the plot to show. If not specified,
        it plots the standard way.

    Returns
    -------
    None.

    """
    
    # Define range
    x_min = sys.float_info.max
    y_min = sys.float_info.max
    
    x_max = (-1) * sys.float_info.max
    y_max = (-1) * sys.float_info.max
    
    for p in points:
        if p[0] < x_min:
            x_min = p[0]
        if p[0] > x_max:
            x_max = p[0]
        
        if p[1] < y_min:
            y_min = p[1]
        if p[1] > y_max:
            y_max = p[1]
            
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    res = 1000
    if x_range < 1e-6:
        x_range = 1e-6
    x_delta = x_range / res
    
    if y_range < 1e-6:
        y_range = 1e-6
    y_delta = y_range / res
    
    x_max += x_range/3
    x_min -= x_range/3
     
    y_max += y_range/3
    y_min -= y_range/3
    
    # Create grid and function matrix
    x = np.arange(x_min, x_max, x_delta)
    y = np.arange(y_min, y_max, y_delta)
    X,Y = np.meshgrid(x,y)
    Z = function([X,Y])
   
    # Plot contours
    if levels == None:
        levels = 49
    if axes != None:
        axes.contour(X, Y, Z, levels, colors=['black'], alpha= 0.8, linewidths = 1)
    else:
        plt.cla()
        plt.contour(X, Y, Z, levels, colors=['black'], alpha= 0.8, linewidths = 1)
    
    # Plot Newton iteration
    p_x = []
    p_y = []
    for p in points:
        p_x.append(p[0])
        p_y.append(p[1])
    if axes != None:
        axes.plot(p_x, p_y, marker="o", color="red", linestyle="dotted", linewidth = 2, markersize=4)
    else:
        plt.plot(p_x, p_y, marker="o", color="red", linestyle="dotted", linewidth = 2, markersize=4)

    
if __name__ == "__main__":
    # Rosenbrock f(x,y) = (1-x)^2 + 100(y-x^2)^2
    f = lambda x: (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
    points = [np.array([-1/2,-1]), np.array([1,2]), np.array([1.5,3])]
    levels = [1, 3.8, 14, 56, 215, 825]
    fig, (ax1, ax2) = plt.subplots(1, 2)
    plot_newton_2d(f, points, levels, ax1)
    plot_newton_2d(f, [np.array([1,2]), np.array([3,4])], levels, ax2)