# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 11:45:47 2020

@author: joelw
"""

import minimization_solver as MS
import minimization_problem as MP
import numpy as np
from hessian_inv import finite_differences, brute_inv_hessian
import matplotlib.pyplot as plt
import chebyquad_problem as cheby
import plotter

if __name__ == '__main__':
    
    
    # Define problem
    f = lambda x: cheby.chebyquad(x)
    g = lambda x: cheby.gradchebyquad(x)
    p = np.linspace(0,1,11)
    problem = MP.minimization_problem(f, p,g)
    
    # Solve minimization
    solver = MS.minimization_solver(problem)
    
    solver.parameter_update(sensitivity = 0.00001,
                             hessian_approximation_method   =  "bfgs",
                             line_search_method             =  "exact",
                             line_search_conditions         =  "WP")
    track = ['xk', 'inv_hessian']
    x,t = solver.solve(track)
    e_h =[]
    dif = []
    # Calculate hessian for each xk
    for i in range(len(t['xk'])):
        exact_hessian =   brute_inv_hessian(f, t['xk'][i])
        e_h.append( exact_hessian)
        dif.append(np.linalg.norm(t['inv_hessian'][i] - exact_hessian ))

    # Plor the difference 
    plt.plot(dif)
    plt.show()
    
    








