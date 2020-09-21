import numpy as np

class minimization_solver():
    
    def __init__(self, minimization_problem, rho = 0.1, sigma = 0.7, tau = 0.1, chi = 9):
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
        self.minimization_problem = minimization_problem
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
    
    # For changing attributes
    def __call__(self, rho = None, sigma = None, tau = None, chi = None):
        """

        Parameters
        ----------
        rho : float, optional
            Method parameter. The default is None.
        sigma : float, optional
            Method parameter. The default is None.
        tau : float, optional
            Method parameter. The default is None.
        chi : float, optional
            Method parameter. The default is None.

        Returns
        -------
        None.

        """
        if rho != None:
            self.rho = rho
        if sigma != None:
            self.sigma = sigma
        if tau != None:
            self.tau = tau
        if chi != None:
            self.chi = chi
        
    def solve(self, hessian_approximation_method = "good_broyden", line_search_method = "exact", line_search_conditions = "goldstein"):
        """

        Parameters
        ----------
        hessian_approximation_method : str, optional
            Method used for approximation of the Hessian matrix
            and its inverse. The default is "good_broyden".
        line_search_method : str, optional
            Method used to decide the step size to take in the Newton direction. 
            The default is "exact".
        line_search_conditions : str, optional
            Conditions for line search in case "inexact" is chosen in the 
            previous parameter. The default is "goldstein".

        Returns
        -------
        None.

        """
        # General Newton method