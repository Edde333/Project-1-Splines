from  scipy import *
from  pylab import *
import numpy as np

def getBaseFunc(u_real, i, k = 3):
    # Solves fact that extra knots seems to be needed but are not
    u = u_real.copy()
    for j in range(k-2):
        u = np.append(u,1)
    return getBaseFuncRec(u,i,k)

def getBaseFuncRec(u, i, k):
    # Recursive algorithm to find base function
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

# Returns s(x)
def getCubicSpline(x, u, d):
    # Checks for correct x's
    if (x > u[-3] + x < u[2]).any():
        raise Exception("Wrong input")
        
    # Finds and saves all base functions
    base_functions = np.array([])
    for i in range(len(u)-2):
        base_functions = np.append(base_functions, getBaseFunc(u, i)
        
    # Finds hot intervals
    I = np.array([])
    for i in range(len(x)):
        for j in range(len(u)):
            if (x[i] >= u[j]) and (x[i] < u[j+1]):
                I = np.append(I, j)
                
    # for i in range(len(u)):
    #     if (x >= u[i]) * (x < u[i+1]):
    #         I = something
    s_x
    for i in range(len(x)):
        control_points_x = np.array([])
        control_points_y = np.array([])
        for j in range(4):
            control_points_x = np.append(control_points_x, d[(I[i]-2) + j][0])
            control_points_y = np.append(control_points_y, d[(I[i]-2) + j][1])
        hot_base_functions = base_functions[I[i]-1 : I[i]+3]
        s_x = control_points_x @ hot_base_functions
        s_y = control_points_y @ hot_base_functions
    
    
    
    return zip(s_x, s_y)