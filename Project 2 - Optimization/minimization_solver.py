
"""
Functions that need to be implemented either in this module or another one.

"""

import numpy as np
from hessian_inv import get_inverse_hessian

def line_search(x,g,iG,
                f,
                line_search_method = "None",
                line_search_conditions = "Goldstein",
                rho = 0.1,
                sigma = 0.7,
                tau = 0.1,
                chi = 9.):
    
    return 1.
    
"""
def get_inverse_Hessian(x,g,iG,
                    f,
                    hessian_approximation_method = "Brute force"):
    
    
    return brute_invHessian(f, x,)
"""    
          



def brute_force_gradient(f, input_dim, h = 0.001):
    """
        Takes in a pyfunction and maybe something 
        of what kind of object that function acts upon
        and then returns the gradient of that function
        as a pyfunction.
        output shape will be the same as the input_shape
    """
    gradient = np.empty(input_dim, dtype = object)
    H = h*np.identity(input_dim)
    for index in range(input_dim):
        gradient[index] = lambda x: (f(x+H[index]) - f(x-H[index]))/(2*h)

    return lambda x: np.array([gradient[i](x) for i in range(input_dim)], dtype = object)

def brute_gradient(f, x, h = 0.0001):
    input_dim = x.size
    gradient = np.empty(input_dim, dtype = float)
    I = np.identity(input_dim)
    for index in range(input_dim):
        gradient[index] = (f(x+h*I[index]) - f(x-h*I[index]))/(2*h)  
    return gradient


def brute_hessian(f, x, h = 0.000001):
    input_dim = x.size
    hessian = np.empty( ( input_dim, input_dim))
    I =np.zeros(input_dim)
    J = np.zeros( input_dim)
    for index in np.ndindex(hessian.shape):
        I[index[0]] = h
        J[index[1]] = h
        hessian[index] = (f(x+I+J)+f(x-I-J) - f(x-I+J) - f(x+I-J))/(4*h*h)
        I[index[0]] = 0
        J[index[1]] = 0
    return hessian

def brute_invHessian(f, x, h = 0.000001):
    hessian = brute_hessian(f, x, h = h)
    hessian = 1/2*(hessian+hessian.T)
    return np.linalg.inv(hessian)

    
def exact_line_search(f, x, g, invG):
    alpha = 0
    f_prev = f(x-alpha*invG@g)
    
    
    
    

   
def exact_line_search(f):
    """
        Takes in a python function and returns a linesearch function 
        that takes in an x_k, s_k, funciton and then iterates over alpha
    """
    pass

def inexact_line_search(f, option):
    """
        Takes in a python function and returns a linesearch function 
        that takes in an x_k, s_k, funciton and then iterates over alpha
        it also takes in an option to choose which condition applies
    """
    pass


class simple_Problem:
    
    def __init__(self, function, guess, gradient = None):
        self.function = function
        self.guess = guess
        if gradient is None:
            self.gradient = lambda x: brute_gradient(function, x)
        else :
            self.gradient = gradient
        
            


# ACTUAL CLASS AND STUFF!!!!


class minimization_solver():
    
    
    
    
    def __init__(self,
                 Problem,
                 f = None,
                 guess = None,
                 gradient = None,
                 hessian_approximation_method = "brute force",
                 line_search_method = "None",
                 line_search_conditions = "Goldstein",
                 rho = 0.1,
                 sigma = 0.7,
                 tau = 0.1,
                 chi = 9,
                 sensitivity = 0.000001
                 ):
        
        #Maybe do some checks of the input first
        # pass
        
        # Initiate the simple variables
        self.Problem = Problem
        #self.f = f
        self.x_k = Problem.guess
        
        #self.gradient = gradient
        self.line_search_method = line_search_method
        self.line_search_conditions = line_search_conditions
        self.hessian_approximation_method = hessian_approximation_method
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
        self.sensitivity = sensitivity
        
        #self.invHess = brute_hessian(self.Problem.function, self.x_k)
        #self.grad = brute_gradient(self.Problem.function, self.x_k)
        
        self.update_methods()
        

    def update_methods(self):
        """
            Updates the methods (gradient, line_seach and inverse_hessian)
            according to the current state of the class
        """
        #if self.gradient is None: 
        #    # Shall return a gradient function vf: R^n -> R^n
        #    self.gradient = lambda x: brute_gradient(self.f, x)
            
        # Handle the options inserted as strings 
        self.next_alpha = lambda x,g,iG: line_search(x,g,iG,
                                                     f = self.Problem.function,
                                                     line_search_method = self.line_search_method ,
                                                     line_search_conditions = self.line_search_conditions,
                                                     rho = self.rho,
                                                     sigma = self.sigma,
                                                     tau = self.tau,
                                                     chi = self.chi)
            
        
        self.next_invHessian = lambda xk, xk_i, prev_hessian: get_inverse_hessian(self.Problem, xk, xk_i, prev_hessian,
                                                              hessian_approximation_method =  self.hessian_approximation_method)


    def __call__(self):
        """
            The method runs the algoritm with the parameters previously set.
            

        Returns
        -------
        x_k : np array, dtype = float
            The algoritms solution to the minimazation problem
            
        """
        
        # Initiating variables
        # Create a local variable
        xk = self.x_k.copy() 
        xk_1 = xk.copy() 
        invHess = get_inverse_hessian(self.Problem.function, xk)
        prev_hessian = get_inverse_hessian(self.Problem.function, xk)
        
        # Do first step
        grad = brute_gradient(self.Problem.function, xk)
        alpha  = self.next_alpha(xk, grad, invHess)
        xk = xk - alpha*invHess@grad
        dx = np.linalg.norm(xk - xk_1)
        while ( dx > self.sensitivity):
            # Calculate inverse hessian, gradiant and alpha
            invHess  =   self.next_invHessian(xk,xk_1, prev_hessian)
            grad     =   self.Problem.gradient(xk)
            alpha    =   self.next_alpha(xk, grad, invHess)
            # update parameters
            xk_1     =   xk.copy()
            x_k      =   xk - alpha*invHess@grad
            dx       =   np.linalg.norm(xk - xk_1)
            print(xk)
        
        return xk

    def parameter_update(self,
                 Problem = None,
                 f = None,
                 guess = None,
                 gradient = None,
                 hessian_approximation_method = None,
                 line_search_method = None,
                 line_search_conditions = None,
                 rho = None,
                 sigma = None,
                 tau = None,
                 chi = None,
                 sensitivity = None
                 ):
        
        self.Problem = Problem if Problem is not None else self.Problem
        self.f = f if f is not None else self.f
        self.x_k = guess.copy() if guess is not None else self.x_k
        
        self.gradient = gradient if gradient is not None else self.gradient
        self.line_search_method = line_search_method if line_search_method is not None else self.line_search_method
        self.line_search_conditions = line_search_conditions if line_search_conditions is not None else self.line_search_conditions
        self.hessian_approximation_method = hessian_approximation_method if hessian_approximation_method is not None else self.hessian_approximation_method
        self.rho = rho if rho is not None else self.rho
        self.sigma = sigma if sigma is not None else self.sigma
        self.tau = tau if tau is not None else self.tau
        self.chi = chi if chi is not None else self.chi
        self.sensitivity = sensitivity if sensitivity is not None else self.sensitivity
        
        self.invHess = brute_hessian(self.f, self.x_k)
        self.grad = brute_gradient(self.f, self.x_k)
        
        self.update_methods()

if __name__ == '__main__':
    f = lambda x: x[0]**2+x[1]**2
    guess = np.array([10,10])
    prob = simple_Problem(f, guess)
    ms = minimization_solver(prob)
    
    xp=ms()
    print(xp)

   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
