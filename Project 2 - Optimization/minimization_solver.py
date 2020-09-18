
"""
Functions that need to be implemented either in this module or another one.

"""

import numpy as np

def line_search(x,g,iG,
                f,
                line_search_method = "None",
                line_search_conditions = "Goldstein",
                rho = 0.1,
                sigma = 0.7,
                tau = 0.1,
                chi = 9.):
    
    return 1.
    

def inverse_Hessian(x,g,iG,
                    f,
                    hessian_approximation_method = "Brute force"):
    
    
    return brute_invHessian(f, x,)
    
          



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

class minimization_solver():
    
    
    
    
    def __init__(self,
                 Problem,
                 f = None,
                 guess = None,
                 gradient = None,
                 hessian_approximation_method = "Brute force",
                 line_search_method = "None",
                 line_search_conditions = "Goldstein",
                 rho = 0.1,
                 sigma = 0.7,
                 tau = 0.1,
                 chi = 9,
                 sensitivity = 0.01
                 ):
        
        #Maybe do some checks of the input first
        # pass
        
        # Initiate the simple variables
        self.Problem = Problem
        self.f = f
        self.x_k = guess
        
        self.gradient = gradient
        self.line_search_method = line_search_method
        self.line_search_conditions = line_search_conditions
        self.hessian_approximation_method = hessian_approximation_method
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
        self.sensitivity = sensitivity
        
        self.invHess = brute_hessian(self.f, self.x_k)
        self.grad = brute_gradient(self.f, self.x_k)
        
        self.update_methods()
        

    def update_methods(self):
        """
            Updates the methods (gradient, line_seach and inverse_hessian)
            according to the current state of the class
        """
        if self.gradient is None: 
            # Shall return a gradient function vf: R^n -> R^n
            self.gradient = lambda x: brute_gradient(self.f, x)
            
        # Handle the options inserted as strings 
        self.next_alpha = lambda x,g,iG: line_search(x,g,iG,
                                                     f = self.f,
                                                     line_search_method = self.line_search_method ,
                                                     line_search_conditions = self.line_search_conditions,
                                                     rho = self.rho,
                                                     sigma = self.sigma,
                                                     tau = self.tau,
                                                     chi = self.chi)
            
        
        self.next_invHessian = lambda x,g,iG: inverse_Hessian(x,g,iG,
                                                              self.f,
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
        dx = 2*self.sensitivity # Just set to larger then the sensetivity
        # Create a local variable
        x_k = self.x_k 
        invHess = self.invHess.copy()
        grad = self.grad.copy()
        
        while ( dx > self.sensitivity):
            # Calculate inverse hessian, gradiant and alpha
            invHess  =   self.next_invHessian(x_k,grad, invHess)
            grad     =   self.gradient(x_k)
            alpha    =   self.next_alpha(x_k, grad, invHess)
            # update parameters
            prev_x   =   x_k.copy()
            x_k      =   x_k - alpha*invHess@grad
            dx       =   np.linalg.norm(x_k-prev_x)
        
        return x_k
        

if __name__ == '__main__':
    f = lambda x: x[0]**2+x[1]**2+x[0]*x[1]+12
    guess = np.array([60,10])
    
    ms = minimization_solver("mini", f = f, guess = guess)
    
    xp=ms()
    print(xp)
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   