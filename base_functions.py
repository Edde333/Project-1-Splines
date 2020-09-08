from  scipy import *
from  pylab import * 
import numpy as np


"""

"""


def getBaseFunc(u_real, i, k = 3):
    # Solves fact that extra knots seems to be needed but are not
    u = u_real.copy()
    if i > len(u_real) - 3 or i < 0:
        raise Exepction("Wrong input")
    
    for j in range(k-2):
        u = np.append(u,1)
    return getBaseFuncRec(u,i,k)



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
        # try:
        #     factor1 = lambda x: (x - u[i-1])/(u[i+k-1] - u[i-1])
        # except ZeroDivisionError:
        #     factor1 = lambda x: 0
            
        # try:
        #     factor2 = lambda x: (u[i+k] - x)/(u[i+k] - u[i])
        # except ZeroDivisionError:
        #     factor2 = lambda x: 0
        
        return lambda x: factor1(x) * getBaseFuncRec(u, i, k-1)(x) + factor2(x) * getBaseFuncRec(u, i+1, k-1)(x)




def getCubicSpline(x, u, d):
    """
    Uses the deBoor algorithm to create splines
    
    Parameters
    ----------
    x : array of floats, array of floats
        the point(s) containd in u, for which the alogorithm computes a spline
    u : array of floats, length K, nodes in u_i
    d : coordinates corresponding to nodes u_1 to u_K-1 


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
    for i in range(len(u)-2):
        base_functions = np.append(base_functions, getBaseFunc(u, i))
    
                                   
    # Finds hot intervals
    I = np.array([]) 
    for i in range(len(x)):
        for j in range(len(u)):
            if (x[i] >= u[j]) and (x[i] < u[j+1]):
                I = np.append(I, j)
                
    # for i in range(len(u)):
    #     if (x >= u[i]) * (x < u[i+1]):
    #         I = something
    
    for i in range(len(x)):
        control_points_z = np.array([])
        control_points_y = np.array([])
        for j in range(4):
            control_points_z = np.append(control_points_z, d[(I[i]-2) + j][0])
            control_points_y = np.append(control_points_y, d[(I[i]-2) + j][1])
        hot_base_functions = base_functions[I[i]-1 : I[i]+3]
        s_z = control_points_z @ hot_base_functions
        s_y = control_points_y @ hot_base_functions

    
    
    return zip(s_z, s_y)




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
    
                                   
    return np.sum(basefunctions(x))