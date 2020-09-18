import numpy as np

class minimization_solver():
    
    def approx_hessian_inv_by_finite_diff(gradient, x):
        """
    
    
        Parameters
        ----------
        gradient : N-array of R^N to R functions
            The gradient of the function for which we want to find a minimum
        x : N-Array
            An array containing the coordinates on which we want to
            approximate the Hessian matrix
    
        Raises
        ------
        Exception
            Raises an Exception if the approximated Hessian is not positive
            definite
    
        Returns
        -------
        NxN-array
            Returns the approximated inverse Hessian matrix for the problem at
            the specified point x
    
        """
        # Test f = 3(x1)**2 + 2*x2
        # x = np.array([1,2])
        # gradient = np.array([lambda x: 6 * x[0] + 4 * x[1],
        #                       lambda x: 4 * x[1] + 4 * x[0]])
        
        dim = len(gradient)
        dx = 1e-6
        hessian = np.zeros(dim**2).reshape(dim,dim)
        for i in range(dim): #row
            for j in range(i,dim): #column
                # creates array used for difference calculation
                pc = np.zeros(dim)
                pc[j] = dx
                
                # carries out derivation
                deriv = (gradient[i](x + pc) - gradient[i](x - pc)) / (2 * dx)
                
                # constructs hessian
                hessian[i,j] = deriv
                hessian[j,i] = deriv
        
        # checks that hessian is positive definite
        try:
            np.linalg.cholesky(hessian)
        except np.linalg.LinAlgError:
            raise Exception("Hessian matrix not positive definite")
            
        return np.linalg.inv(hessian)
            
            
