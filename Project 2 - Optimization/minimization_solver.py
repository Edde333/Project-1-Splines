
"""
Functions that need to be implemented either in this module or another one.

"""

def brute_force_gradient(f):
    """
        Takes in a pyfunction and maybe something 
        of what kind of object that function acts upon
        and then returns the gradient of that function
        as a pyfunction.
    """
    # Below shall change
    return 

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
                 gradient = None,
                 hessian_approximation_method = "Brute force",
                 line_search_method = "None",
                 line_search_conditions = "Goldstein",
                 rho = 0.1,
                 sigma = 0.7,
                 tau = 0.1,
                 chi = 9,
                 sensitivity = 0.00000000001
                 ):
        
        #Maybe do some checks of the input first
        # pass
        
        # Initiate the simple variables
        self.Problem
        self.x_k = Problem.guess.copy()
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
        
        
        
        # Handle gradient which might be included
        if gradient is None: 
            # Shall return a gradient function vf: R^n -> R^n
            gradient = brute_force_gradient(Problem.function)
            
        # Handle the options inserted as strings 
        if line_search_method == "None":
            # Something easy...
            get_next_alpha = lambda x: 1.
        elif line_search_method == "Exact":
            # I guess we had to pass on the condition option here?
            # There might be a nicer way around it...
            # Maybe the stepside shall be passed here?
            get_next_alpha = exact_line_search(Problem.function)
        elif line_search_method == "Inexact":
            # I guess we had to pass on the condition option here?
            # There might be a nicer way around it...
            # Maybe the stepside shall be passed here?
            option = "Something from the input... Conditions?"
            get_next_alpha = inexact_line_search(Problem.function, option)
            
        if hessian_approximation_method == "Brute force":
            """
                It works the same here as in the if statments above
            """
            get_next_invHessian = "somefunction"
        elif hessian_approximation_method == "one of 5 options...":
            pass

    def __call__(self):
        dx = 2*self.sensitivity
        x_k = self.x_k
        while ( dx < self.sensitivity):
            invHess  =   self.get_next_invHessian(x_k)
            grad     =   self.gradient(self.x_k)
            alpha    =   self.get_next_alpha(x_k)
            prev_x   =   x_k.copy()
            x_k      =   x_k - alpha*invHess*grad
            dx       =   (x_k-prev_x).norm()
        
        return x_k
            
        
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   
    
   