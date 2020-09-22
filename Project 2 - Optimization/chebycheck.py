
"""
Created on Fri Sep 18 10:50:01 2020

@author: kadde
"""

import chebyquad_problem as cheby
import minimization_problem as mini
import minimization_solver as slv
import numpy as np
import scipy as sp
import scipy.optimize as so

f = lambda x: cheby.chebyquad(x)
g = lambda x: cheby.gradchebyquad(x)

guess = np.linspace(0,1,8)
problem = mini.minimization_problem(f,guess)
solver = slv.minimization_solver(problem,'bfgs', 'inexact','WP')
xmin1 = solver.solve()


xmin= so.fmin_bfgs(f,guess,g)  # should converge after 18 iterations 
print(xmin)
print(xmin1)
print(f(xmin))
print(f(xmin1))
