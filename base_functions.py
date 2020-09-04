from  scipy import *
from  pylab import *

def getBaseFunc(u, i, k = 3):
    if k == 0:
        if u[i-1] == u[i]:
            return lambda x: 0
        else:
            return lambda x: 1 if (x >= u[i-1] and x < u[i]) else 0
    else:
        try:
            factor1 = lambda x: (x - u[i-1])/(u[i+k-1] - u[i-1])
        except ZeroDivisionError:
            factor1 = lambda x: 0
            
        try:
            factor2 = lambda x: (u[i+k] - x)/(u[i+k] - u[i])
        except ZeroDivisionError:
            factor2 = lambda x: 0
        
        return lambda x: factor1(x) * getBaseFunc(u, i, k-1)(x) + factor2(x) * getBaseFunc(u, i+1, k-1)(x)
