import numpy as np
from hessian_inv import get_inverse_hessian, get_accaptable_hessian_approximation_methods
from line_search import line_search
from minimization_problem import minimization_problem
import matplotlib.pyplot as plt
import sys


class minimization_solver():

    def __init__(self, minimization_problem,
                 hessian_approximation_method   =  "finite_differences",
                 line_search_method             =  "exact",
                 line_search_conditions         =  "goldstein",
                 rho    =  0.1,
                 sigma  =  0.7,
                 tau    =  0.1,
                 chi    =  9.,
                 sensitivity  =  0.1):
        """
        Initializes a solver.

        Parameters
        ----------
        minimization_problem : minimization_problem
            Problem to be solved by this solver.
        rho : float, optional
            Method parameter. The default is 0.1.
        sigma : float, optional
            Method parameter. The default is 0.7.
        tau : float, optional
            Method parameter. The default is 0.1.
        chi : float, optional
            Method parameter. The default is 9.

        Returns
        -------
        None.

        """
        self.minimization_problem           = minimization_problem
        self.hessian_approximation_method   = hessian_approximation_method
        self.line_search_method             = line_search_method
        self.line_search_conditions         = line_search_conditions
        self.rho                            = rho
        self.sigma                          = sigma
        self.tau                            = tau
        self.chi                            = chi
        self.sensitivity                    = sensitivity
        
        self._check_parameter()

    # For changing attributes
    def parameter_update(self,
                 minimization_problem          = None,
                 hessian_approximation_method  = None,
                 line_search_method            = None,
                 line_search_conditions        = None,
                 rho                           = None,
                 sigma                         = None,
                 tau                           = None,
                 chi                           = None,
                 sensitivity                   = None):
        """
        Parameters
        ----------
                 minimization_problem          : minimization_problem object. 
                 hessian_approximation_method  : string, sets the method of how the hessian is calculated.
                 line_search_method            : string, sets which method is used as line search.
                 line_search_conditions        : string, sets which condition is used in line search.
                 rho                           : float, parameter used in linesearch.
                 sigma                         : float, parameter used in linesearch.
                 tau                           : float, parameter used in linesearch.
                 chi                           : float, parameter used in linesearch.
                 sensitivity                   : float, the maximum distance that x is updated before algorithm termination.
        
        Updates parameters in the solver object which dictates which problem to be solved and how it is to be solved.
        If no value is given, the parameter is not changed. If none-accaptable input is given, an exception is raised.
        Returns
        -------
        None.

        """
        
        self.minimization_problem            = minimization_problem          if minimization_problem is not None else self.minimization_problem
        self.line_search_method              = line_search_method            if line_search_method is not None else self.line_search_method
        self.line_search_conditions          = line_search_conditions        if line_search_conditions is not None else self.line_search_conditions
        self.hessian_approximation_method    = hessian_approximation_method  if hessian_approximation_method is not None else self.hessian_approximation_method
        self.rho                             = rho                           if rho is not None else self.rho
        self.sigma                           = sigma                         if sigma is not None else self.sigma
        self.tau                             = tau                           if tau is not None else self.tau
        self.chi                             = chi                           if chi is not None else self.chi
        self.sensitivity                     = sensitivity                   if sensitivity is not None else self.sensitivity
        
        self._check_parameter()

    def solve(self, tracking_names = []):
        """
        Parameters
        ----------
        tracker_names : list of strings.
            Contains the names of the variables you want to check.
            accaptable values: ['xk', 'xk_1', 'inv_hessian', 'grad', 'alpha', 'dx']

        Returns
        -------
        xk: n-array
            The algoritms choice as a minimum point.
        tracker: dict: str: list
            Contains the list of the progression of each
        """
        # Check if input is acceptable
        acceptable_names = ['xk', 'xk_1', 'inv_hessian', 'grad', 'alpha', 'dx']
        if not set(tracking_names).issubset(acceptable_names):
            raise Exception("tracking_names must be a subset of {}".format(acceptable_names))
        
        
        # Initiating variables
        tracker       = {name : [] for name in tracking_names}
        xk_1          = None
        xk            = self.minimization_problem.guess.copy()
        dx            = None
        grad          = self.minimization_problem.gradient(xk)
        inv_hessian   = get_inverse_hessian(self.minimization_problem, 
                                            xk,
                                            None,
                                            None,
                                            hessian_approximation_method = "finite_differences")
        
        alpha         = line_search(self.minimization_problem.function,
                                    self.minimization_problem.gradient,
                                    xk,
                                    -inv_hessian@grad,
                                    line_search_method = self.line_search_method,
                                    line_search_condition = self.line_search_conditions, 
                                    a0 = 1, rho = 0.1, sigma = 0.7, tau = 0.1, chi = 9)
        
        # Save initial values
        for name,l in tracker.items(): l.append(locals()[name])
        
        # Continue with next step according to object settings
        while (True):
            # Update parameters
            xk_1     =   xk
            xk       =   xk - alpha*inv_hessian@grad
            dx       =   np.linalg.norm(xk - xk_1)
            # Calculate the proerties at new xk
            grad           =   self.minimization_problem.gradient(xk)
            inv_hessian    =   get_inverse_hessian(self.minimization_problem,
                                                 xk, 
                                                 xk_1,
                                                 inv_hessian,
                                                 self.hessian_approximation_method)

            alpha          =   line_search(self.minimization_problem.function,
                                           self.minimization_problem.gradient,
                                           xk,
                                           -inv_hessian@grad,
                                           line_search_method = self.line_search_method,
                                           line_search_condition = self.line_search_conditions, 
                                           a0 = 1, rho = self.rho, sigma = self.sigma, tau = self.tau, chi = self.chi)

            # Save values
            for name,l in tracker.items(): l.append(locals()[name].copy())
            
            # Check if condition is reached
            #if dx <   self.sensitivity: break
            if np.linalg.norm(grad) < self.sensitivity: break
            #if abs(self.minimization_problem.function(xk)- self.minimization_problem.function(xk_1)) < self.sensitivity: break
        
        return xk, tracker
    
    def _check_parameter(self):
        """
        Raises
        ------
        
            Exception: If the set parameters are none-accaptable

        Returns
        -------
        None.

        """
        accaptable_hessian_approximation_method   =  get_accaptable_hessian_approximation_methods()  
        accaptable_line_search_method             =  ['exact', 'inexact']
        accaptable_line_search_conditions         =  ["goldstein", 'WP']
        
        # Check if accaptable
        if self.hessian_approximation_method not in accaptable_hessian_approximation_method:
            raise  Exception("hessian approximation method must be one in:\n{}\nGiven value: {}".format(accaptable_hessian_approximation_method,
                                                                                     self.hessian_approximation_method))
        if self.line_search_method not in accaptable_line_search_method:
            raise  Exception("line search method must be one in:\n{}\nGiven value: {}".format(accaptable_line_search_method,
                                                                                              self.line_search_method))
        if self.line_search_conditions not in accaptable_line_search_conditions:
            raise  Exception("line search conditions must be one in {}\nGiven value: {}".format(accaptable_line_search_conditions,
                                                                                                 self.line_search_conditions))
        """   
        rho    =  0.1,
        sigma  =  0.7,
        tau    =  0.1,
        chi    =  9.,
        sensitivity  =  0.01
        """
       
    def plot_newton_2d(self, points, levels = None):
        
        """
        Parameters
        ----------
        points : (Nx2)-array
            The recursive points that you wish to plot.
        levels : array-like, optional
            An array with the (increasing) levels that the contour plot
            is drawing. If not specified 50 equally spaced (value-wise)
            contours are drawn. Default is None.
    
        Returns
        -------
        None.
    
        """
        function = self.minimization_problem.function
        # Define range
        x_min = sys.float_info.max
        y_min = sys.float_info.max
        
        x_max = sys.float_info.min
        y_max = sys.float_info.min
        
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
        x_max += x_range/3
        x_min -= x_range/3
        
        y_range = y_max - y_min
        y_max += y_range/3
        y_min -= y_range/3
        
        # Create grid
        delta = 0.025
        x = np.arange(x_min, x_max, delta)
        y = np.arange(y_min, y_max, delta)
        X,Y = np.meshgrid(x,y)
        Z = function([X,Y])
       
        
        if levels == None:
            levels = 49
        plt.contour(X, Y, Z, levels, colors=['black'], alpha= 0.8, linewidths = 1)
        
        p_x = []
        p_y = []
        for p in points:
            p_x.append(p[0])
            p_y.append(p[1])
            
        plt.plot(p_x, p_y, marker="o", color="red", linestyle="dotted", linewidth = 2, markersize=4)
        plt.show()
                
   
if __name__ == '__main__':
    import math
    f = lambda x: x[0]**2+x[1]**2
    guess = np.array([3.8,8])
    problem = minimization_problem(f,guess)
    solver = minimization_solver(problem)
    solver.parameter_update(hessian_approximation_method="good_broyden",
                            line_search_method = "inexact",)
    names = ['xk', 'inv_hessian', 'xk_1']
    x,tracker = solver.solve(names)
    solver.plot_newton_2d(tracker['xk'])