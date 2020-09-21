import numpy as np
from hessian_inv import get_inverse_hessian
#from line_search import line_search
from minimization_problem import minimization_problem

class minimization_solver():

    def __init__(self, minimization_problem, hessian_approximation_method = "brute_inverse_hessian",
                                             line_search_method = "exact",
                                             line_search_conditions = "goldstein",
                                             rho = 0.1,
                                             sigma = 0.7,
                                             tau = 0.1,
                                             chi = 9,
                                             sensitivity =0.01):
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
        self.hessian_approximation_method = hessian_approximation_method
        self.line_search_method = line_search_method
        self.line_search_conditions = line_search_conditions
        self.rho = rho
        self.sigma = sigma
        self.tau = tau
        self.chi = chi
        self.sensitivity = sensitivity

    # For changing attributes
    def parameter_update(self,
                 minimization_problem = None,
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

        self.minimization_problem            = minimization_problem if minimization_problem is not None else self.minimization_problem
        self.line_search_method              = line_search_method if line_search_method is not None else self.line_search_method
        self.line_search_conditions          = line_search_conditions if line_search_conditions is not None else self.line_search_conditions
        self.hessian_approximation_method    = hessian_approximation_method if hessian_approximation_method is not None else self.hessian_approximation_method
        self.rho                             = rho if rho is not None else self.rho
        self.sigma                           = sigma if sigma is not None else self.sigma
        self.tau                             = tau if tau is not None else self.tau
        self.chi                             = chi if chi is not None else self.chi
        self.sensitivity                     = sensitivity if sensitivity is not None else self.sensitivity

    def solve(self):
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
              # Initiating variables
        # Create a local variable
        xk = self.minimization_problem.guess.copy()
        xk_1 = xk.copy()
        inv_hessian = get_inverse_hessian(self.minimization_problem, xk, None, None, hessian_approximation_method = "brute_inverse_hessian")

        # Do first step
        grad = self.minimization_problem.gradient(xk)
        #alpha  = line_search(...)
        alpha  = 1.
        xk = xk - alpha*inv_hessian@grad
        dx = np.linalg.norm(xk - xk_1)
        while ( dx > self.sensitivity):
            # Calculate inverse hessian, gradiant and alpha
            inv_hessian  =   get_inverse_hessian(self.minimization_problem, xk, xk_1, inv_hessian, self.hessian_approximation_method)
            grad     =   self.minimization_problem.gradient(xk)
            #alpha  = line_search(...)
            alpha    =   1.
            # update parameters
            xk_1     =   xk
            xk       =   xk - alpha*inv_hessian@grad
            dx       =   np.linalg.norm(xk - xk_1)

        return xk


if __name__ == '__main__':
    f = lambda x: x[0]**2+5*x[1]**2+23+10*x[1]
    guess = np.array([12,20])
    problem = minimization_problem(f,guess)
    solver = minimization_solver(problem)
    x = solver.solve()
    print(x)