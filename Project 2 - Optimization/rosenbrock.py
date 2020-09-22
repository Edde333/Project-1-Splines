from minimization_problem import minimization_problem
from minimization_solver import minimization_solver
from plotter import plot_newton_2d
import numpy as np
from matplotlib import pyplot as plt

def check_rosenbrock(guess,
                     hessian_approximation_method = "finite_differences"):
    
    # Task 5 - Rosenbrock with classical Newton and exact line search
    rosenbrock = lambda x: (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
    rosenbrock_grad = lambda x: np.array([2*(x[0]-1)-400*x[0]*(x[1]-x[0]**2),
                                          200*(x[1]-x[0]**2)])
    mp = minimization_problem(rosenbrock, guess, rosenbrock_grad)
    
    ms = minimization_solver(mp, hessian_approximation_method, "exact",
                             #sensitivity = 0.001,
                             #sigma = 0.1,
                             #rho = 0.01
                             )
    res, points = ms.solve(['xk'])
    points = points.get('xk')
    levels = [1, 3.8, 14, 56, 215, 825]
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.set_title("Exact line search")
    plot_newton_2d(ax1, rosenbrock, points, levels)
    print("Exact line search, minimum at: ", res)
    
    # Task 7 - Rosenbrock with classical Newton and inexakt line search
    
    ms.parameter_update(line_search_method="inexact")
    res, points = ms.solve(['xk'])
    points = points.get('xk')
    ax2.set_title("Inexact line search")
    plot_newton_2d(ax2, rosenbrock, points, levels)
    print("Inexact line serach, minimum at: ", res)
    
if __name__ == "__main__":
    guess = np.array([0, -0.7])
    check_rosenbrock(guess, "finite_differences")