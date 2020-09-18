import numpy as np

class minimization_solver():
    
    def find_initial_hessian_inv(min_problem):
        gradient = min_problem.gradient.copy()
        guess = min_problem.guess.copy()
        
        # Test f = 3(x1)**2 + 2*x2
        # guess = np.array([1,2])
        # gradient = np.array([lambda x: 6 * x[0] + 4 * x[1],
        #                      lambda x: 2 + 4 * x[0]])
        
        dim = len(gradient)
        dx = 1e-4 
        hessian = np.zeros(dim**2).reshape(dim,dim)
        for i in range(dim): #row
            for j in range(i,dim): #column
                # creates array used for difference calculation
                pc = np.zeros(dim)
                pc[j] = dx
                
                # carries out derivation
                deriv = (gradient[i](guess + pc) - gradient[i](guess - pc)) / (2 * dx)
                
                # constructs hessian
                hessian[i,j] = deriv
                hessian[j,i] = deriv
        
        return np.linalg.inv(hessian)
            
            
