import numpy as np
import scipy as sp
from scipy import optimize

"""
Method to find an appropriate step size alpha, using the provided line search method



Parameters
----------
f   :   function 

g   :   function, the gradient of f

xk  :   np array, dtype = float, the active coordinates

sk  :   np array, dtype = float, direction for minimization given coordinates xk

line_search_method   : string, specifies which line search method to be used, 
                        accepts "exact" and "inexact", default value "exact"
                        
line_searc_condition : string, specifies which conditions to be used for 
                        inexact line search, accepts "goldstein" and "WP", 
                        default value "goldstein"

a0  :   float, current step size, default value 1             

rho :   float, constant acceptable values between 0 and 0.5
        default value 0.1

sigma : float, constant acceptable values between 0 and 1, sigma >= rho
        default value 0.7

tau :   float, constant
        default value 0.1

chi :   float, constant
        default value 9



Raises
------
    ValueError
        Raises a ValueError if rho or sigma are unacceptable input values
    Invalid line search method
    

Returns
-------
a0 : float
        Step size 
"""


def line_search(f, g, xk, sk, line_search_method = "exact", line_search_condition = "goldstein", 
                a0 = 1, rho = 0.1, sigma = 0.7, tau = 0.1, chi = 9):
    
    
    
    
    if line_search_method == 'none':
        return 1.
    
    
    
    elif line_search_method == "exact":
        fa = lambda a: f(np.add(xk, a * sk))    
        alpha = sp.optimize.fmin(fa, 0, disp=False)

        return alpha
            
              


    elif line_search_method == "inexact":
        
        if (rho < 0 or rho > 0.5):
            raise ValueError("Input rho out of bounds, acceptable values 0 < rho < 0.5")
        if (sigma < 0 or sigma > 1 or sigma < rho):
            raise ValueError("Input rho out of bounds, acceptable values 0 < sigma < 1 with sigma >= rho")
        
        aL = 0
        aU = 10**99
        a0k = (np.add(xk, a0 * sk))
        aLk = (np.add(xk, aL * sk))

        fa0 = f(a0k)
        faL = f(aLk)
        ga0 = g(a0k).T @ sk
        gaL = g(aLk).T @ sk        
        
        [LC, RC] = get_conditions(line_search_condition, a0, aL, aU, rho, sigma, tau, chi, fa0, faL, ga0, gaL)
        
        while not (LC and RC):
           
            if not LC:
                da0 = extrapolate(a0, aL, ga0, gaL)
                da0 = max(da0, tau * (a0 - aL))
                da0 = min(da0, chi * (a0 - aL))
                aL = a0
                a0 = a0 + da0
                
            elif LC :
                aU = min(a0, aU)
                abar = interpolate(a0, aL, fa0, faL, gaL)
                abar = max(abar, aL + tau * (aU - aL))
                abar = min(a0, aU - tau * (aU - aL))
                a0 = abar
            
            a0k = (np.add(xk, a0 * sk))
            aLk = (np.add(xk, aL * sk))

            fa0 = f(a0k)
            faL = f(aLk)
            ga0 = g(a0k).T @ sk
            gaL = g(aLk).T @ sk
                        
            [LC, RC] = get_conditions(line_search_condition, a0, aL, aU, rho, sigma, tau, chi, fa0, faL, ga0, gaL)
        
        return a0
    
    else:
        raise Exception("Invalid line search method")
    
 
    
 
    
"""
Method to check if the stepsize a0 fulfills the inexact line search conditions


Parameters
----------
                        
line_searc_condition : string, specifies which conditions to be used for 
                        inexact line search, accepts "goldstein" and "WP"
a0  :   float, current step size  
aL  :   float, lower end of step size interval 
aU  :   float, upper end of step size interval        
rho :   float, constant 
sigma : float, constant
tau :   float, constant
chi :   float, constant
fa0 :   float, function value given step size a0
faL :   float, function value given step size aL
ga0 :   float, gradient value given step size a0
gaL :   float, gradient value given step size aL



Raises
------
Invalid line search condition


Returns
-------
[LC, RC] : list of booleans
        LC is the left condition 
        RC is the right condition 
"""    
def get_conditions(line_search_condition, a0, aL, aU, rho, sigma, tau, chi, fa0, faL, ga0, gaL):
    
    
    if line_search_condition == "goldstein":

        LC = fa0 >= (faL + (1 - rho) * (a0-aL) * gaL)
        #print(" LC = ", fa0 , ">=",  (faL + (1 - rho) * (a0-aL) * gaL), LC)
        RC = fa0 <= (faL + rho * (a0 - aL) * gaL.T)
        #print("RC = ", fa0," <= ", (faL + rho * (a0 - aL) * gaL), RC)
        return [LC, RC]
        
    elif line_search_condition == "WP":
        LC = (ga0  >= sigma * gaL) 
        #print(" LC = " , ga0, " >= ", sigma * gaL, LC)
        RC = (fa0 <= (faL + rho * (a0 - aL) * gaL.T))
        #print ("RC = ", fa0 ," <= ", (faL + rho * (a0 - aL) * gaL), RC)
        return [LC, RC]        
    
    else: 
        raise Exception("Invalid line search condition")
    




"""
Method to calculate the extrapolated da0, used in inexact line search


Parameters
----------
a0  :   float, current step size  
aL  :   float, lower end of step size interval  
ga0 :   float, gradient value given step size a0
gaL :   float, gradient value given step size aL 

Returns 
-------
da0 : float 
"""
    
def extrapolate(a0, aL, ga0, gaL):
  
    return (a0 - aL) * ga0 / (gaL - ga0)

    


"""
Method to calculate the interpolated abar, used in inexact line search

Parameters
----------
a0  :   float, current step size  
aL  :   float, lower end of step size interval  
fa0 :   float, function value given step size a0
faL :   float, function value given step size aL
gaL :   float, gradient value given step size aL 

Returns 
-------
abar : float 
""" 
def interpolate(a0, aL, fa0, faL, gaL):

    return (a0 - aL)**2 * gaL / (2 * (faL - fa0 + (a0 - aL) * gaL))



"""
Method to return acceptable line search methods, and acceptable 
line search conditions for inexact line search

Returns
-------
acceptable_list_methods, acceptable_list_conditions
"""

def get_acceptable_line_search():
    acceptable_list_methods = [ "exact",
            "inexact",
            "none",
            ]
    acceptable_list_conditions = ["goldstein", "WP"]
    return acceptable_list_methods, acceptable_list_conditions
