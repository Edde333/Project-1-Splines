import numpy as np

class minimization_solver():
    
    def __init__(self, minimization_problem, rho = 0.1, sigma = 0.7, tau = 0.1, chi = 9):
        self.minimization_problem = minimization_problem
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
    
    # For changing attributes
    def __call__(self, rho, sigma, tau, chi):
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
        
    def solve(self, hessian_approximation_method = "good_broyden", line_search_method = "exact", line_search_conditions = "goldstein"):
        
        # General Newton method