import numpy as np

def get_inverse_hessian(minimization_problem, xk, xk_1, prev_hessian, hessian_approximation_method = "good_broyden"):
    """
    

    Parameters
    ----------
    minimization_problem : minimization_problem
        An object containing information on the problem which we want to find
        a minimum
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method
    prev_hessian: NxN-array
        The previous inverse hessian matrix. 
    hessian_approximation_method : String, optional        
        A string specifying the type of hessian inversion method which will be
        used. Valid values of this parameters are:
            * 'good_broyden'
            * 'bad_broyden'
            * 'symmetric_broyden'
            * 'dfp'
            * 'bfgs'
            
        The default is "good_broyden".

    Raises
    ------
    ValueError
        Raises a ValueError if an invalid method is specified in 
        hessian_approximation_method.

    Returns
    -------
    NxN-array
        Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    """
    function = minimization_problem.function
    gradient = minimization_problem.gradient
    
    if hessian_approximation_method == "good_broyden":
        return good_broyden(function, gradient, xk, xk_1, prev_hessian)
    elif hessian_approximation_method == "bad_broyden":
        return bad_broyden(function, gradient, xk, xk_1, prev_hessian)
    elif hessian_approximation_method == "symmetric_broyden":
        return symmetric_broyden(function, gradient, xk, xk_1, prev_hessian)
    elif hessian_approximation_method == "dfp":
         return DFP(function, gradient, xk, xk_1, prev_hessian)
    elif hessian_approximation_method == "bfgs":
         return BFGS(function, gradient, xk, xk_1, prev_hessian)
    else:
        output = "There is no method called: " + hessian_approximation_method
        raise ValueError(output)
            
            
            

def good_broyden(function, gradient, xk, xk_1, prev_hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the good Broyden-method

    Parameters
    ----------
    function : Lamda function
        The function that should be minimized.
    gradient : Lambda function
        The gradient of the function that should be minimized.
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method
    prev_hessian : NxN-array
        The previous hessian matrix.    
        
    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    '''
    deltak = xk - xk_1
    gammak = gradient(xk) - gradient(xk_1)

    new_hess = prev_hessian + (np.outer((deltak - prev_hessian@gammak),(prev_hessian@deltak)))/(np.inner((prev_hessian@deltak),gammak))
    return new_hess


def bad_broyden(function, gradient, xk, xk_1, prev_hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the bad Broyden-method

    Parameters
    ----------
    function : Lamda function
        The function that should be minimized.
    gradient : Lambda function
        The gradient of the function that should be minimized.
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method
    prev_hessian : NxN-array
        The previous hessian matrix.    

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    '''
    deltak = xk - xk_1
    gammak = gradient(xk) - gradient(xk_1)
    
    new_hess = prev_hessian + np.outer(((deltak - prev_hessian@gammak)/(np.inner(gammak,gammak))),gammak)
    return new_hess

def symmetric_broyden(function, gradient, xk, xk_1, prev_hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the symmetric Broyden-method

    Parameters
    ----------
    function : Lamda function
        The function that should be minimized.
    gradient : Lambda function
        The gradient of the function that should be minimized.
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method
    prev_hessian : NxN-array
        The previous hessian matrix.    
        
    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    '''
    pass

def DFP(function, gradient, xk, xk_1, prev_hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the DFP-method

    Parameters
    ----------
    function : Lamda function
        The function that should be minimized.
    gradient : Lambda function
        The gradient of the function that should be minimized.
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method
    prev_hessian : NxN-array
        The previous hessian matrix.    
        
    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the DFP-method.

    '''
    pass

def BFGS(function, gradient, xk, xk_1, prev_hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the BFGS-method

    Parameters
    ----------
    function : Lamda function
        The function that should be minimized.
    gradient : Lambda function
        The gradient of the function that should be minimized.
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method
    prev_hessian : NxN-array
        The previous hessian matrix.    

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the BFGS-method

    '''
    pass
        

if __name__ == "__main__":
    function = lambda x: x[0]**3 + x[1]**3
    gradient = lambda x: np.array([6*x[0],6*x[1]])
    xk_1 = np.array([1, 1])
    xk = np.array([2, 2])
    prev_hessian = np.array([[6, 0],[0,6]])
    y = function(xk)
    new_hessian = good_broyden(function,gradient, xk, xk_1, prev_hessian)
    print(new_hessian)
