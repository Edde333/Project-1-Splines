from minimization_problem import minimization_problem
from minimization_solver import minimization_solver
from plotter import plot_newton_2d
import numpy as np

if __name__ == "__main__":
    # Task 5 - Rosenbrock with classical Newton and exact line search
    rosenbrock = lambda x: (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
    rosenbrock_grad = lambda x: np.array([2*(x[0]-1)-400*x[0]*(x[1]-x[0]**2),
                                          200*(x[1]-x[0]**2)])
    mp = minimization_problem(rosenbrock, [0,-0.7], rosenbrock_grad)
    ms = minimization_solver(mp, "finite_differences", "exact")
    res, points = ms.solve(['xk'])
    points = points.get('xk')
    levels = [1, 3.8, 14, 56, 215, 825]
    plot_newton_2d(rosenbrock, points, levels)
    
    # 