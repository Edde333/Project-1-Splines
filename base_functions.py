from  scipy import *
from  pylab import * 
import numpy as np


"""

"""


def getBaseFunc(u_real, i):
    # Solves fact that extra knots seems to be needed but are not
    u = u_real.copy()
    u = np.append(u,1)

    if i > len(u_real) - 3 or i < 0:
        raise Exception("Wrong input")
    
    # Handles the last base function separately (Defines it as 1 in endpoint)
    if i == len(u_real) - 3:
        return lambda x: (x == u_real[-3]) + getBaseFuncRec(u,i)(x)
    else:
        return getBaseFuncRec(u,i)

# Recursive algorithm to find base function
def getBaseFuncRec(u, i, k = 3):
    """
    Recursive method that finds the base function corresponding with 
    a certain index in the u-vector


    Parameters 
    --------
    u : array of floats, nodes in u_i
    i : integer, the index of the relevant node u_i
    k : integer, the degree of the desired polynomial, is decreased with every iteration 


    Returns 
    -------
    Function (x)
        The base function as a polynomial of degree k
    
    """  
    if k == 0:
        if u[i-1] == u[i]:
            return lambda x: 0
        else:
            return lambda x: (x >= u[i-1]) * (x < u[i])
    else:
        if (u[i+k-1] - u[i-1]) == 0:
            factor1 = lambda x: 0
        else:
            factor1 = lambda x: (x - u[i-1])/(u[i+k-1] - u[i-1])
            
        if (u[i+k] - u[i]) == 0:
            factor2 = lambda x: 0
        else:
            factor2 = lambda x: (u[i+k] - x)/(u[i+k] - u[i])
        
        return lambda x: factor1(x) * getBaseFuncRec(u, i, k-1)(x) + factor2(x) * getBaseFuncRec(u, i+1, k-1)(x)


def getCubicSpline(x, u, d):
    """
    Uses the deBoor algorithm to create splines
    
    Parameters
    ----------
    x : array of floats, array of floats
        the point(s) containd in u, for which the alogorithm computes a spline
    u : array of floats, length K, nodes in u_i
    d : coordinates corresponding to nodes u_0 to u_K-2 


    Raises
    ------
    Exception "Wrong input"
        if any index x lays outside the bounds of u

    Returns
    -------
    (s(z), s(y))
        Tuple with coordinates for point(s) x
    """
    
    I = np.array([]) 
    # Checks if all points are contained within u
    if (x > u[-3]).any() or (x < u[2]).any():
        raise Exception("Wrong input")
    
        
    # Finds and saves all base functions in array base_functions
    base_functions = np.array([])
    for i in range(len(u)-1):
        base_functions = np.append(base_functions, getBaseFunc(u, i))
                                   
    # Finds hot intervals
    I = []
    for i in range(len(x)-1):
        for j in range(len(u)):
            if (x[i] >= u[j]) and (x[i] < u[j+1]):
                I.append(j)
    
    s_z = np.array([])
    s_y = np.array([])
    
    for i in range(len(x)-1):
        control_points_z = np.array([])
        control_points_y = np.array([])
        for j in range(4):
            control_points_z = np.append(control_points_z, d[(I[i]-2) + j][0])
            control_points_y = np.append(control_points_y, d[(I[i]-2) + j][1])
        hot_base_functions = base_functions[I[i]-2 : I[i]+2]
        bf_result = np.array([])
        for bf in hot_base_functions:
            bf_result = np.append(bf_result, bf(x[i]))
        s_z = np.append(s_z, control_points_z @ bf_result)
        s_y = np.append(s_y, control_points_y @ bf_result)
    
    # Last control point is multiplied by 1 in endpoint
    s_z = np.append(s_z, d[-1][0])
    s_y = np.append(s_y, d[-1][1])

    return np.array((s_z, s_y)).T

def baseFuncSum(x, u):
    """
    Returns the sum of all base functions in points x

    Parameters
    ----------
    x : array of floats, points for which the algorithm should work
    u : array of floats, nodes in u_i

    Returns
    -------
    array of floats, sum of all base functions in the given points x

    """
    base_functions = np.array([])
    for i in range(len(u)-2):
        base_functions = np.append(base_functions, getBaseFunc(u, i))
    
    result = np.array([])
    for i in range(len(x)):
        sum = 0
        for bf in base_functions:
            sum += bf(x[i])
        result = np.append(result, sum)
    
    return result