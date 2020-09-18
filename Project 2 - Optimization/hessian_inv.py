

def get_hessian_inv(minimization_problem, xk, xk_1, hessian, hessian_approximation_method = "good_broyden"):
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
    hessian: NxN-array
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
    
    if hessian_approximation_method == "good_broyden":
        return good_broyden(minimization_problem, xk, xk_1, hessian)
    elif hessian_approximation_method == "bad_broyden":
        return bad_broyden(minimization_problem, xk, xk_1, hessian)
    elif hessian_approximation_method == "symmetric_broyden":
        return symmetric_broyden(minimization_problem,xk, xk_1, hessian)
    elif hessian_approximation_method == "dfp":
         return DFP(minimization_problem, xk ,xk_1, hessian)
    elif hessian_approximation_method == "bfgs":
         return BFGS(minimization_problem, xk, xk_1, hessian)
    else:
        output = "There is no method called: " + hessian_approximation_method
        raise ValueError(output)
            
            
            
            

def good_broyden(minimization_problem, xk, xk_1, hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the good Broyden-method

    Parameters
    ----------
    minimization_problem : minimization_problem
        An object containing information on the problem which we want to find
        a minimum
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    '''
    pass

def bad_broyden(minimization_problem, xk, xk_1, hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the bad Broyden-method

    Parameters
    ----------
    minimization_problem : minimization_problem
        An object containing information on the problem which we want to find
        a minimum
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    '''
    pass

def symmetric_broyden(minimization_problem, xk, xk_1, hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the symmetric Broyden-method

    Parameters
    ----------
    minimization_problem : minimization_problem
        An object containing information on the problem which we want to find
        a minimum
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1

    '''
    pass

def DFP(minimization_problem, xk, xk_1, hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the DFP-method

    Parameters
    ----------
    minimization_problem : minimization_problem
        An object containing information on the problem which we want to find
        a minimum
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the DFP-method.

    '''
    pass

def BFGS(minimization_problem, xk, xk_1, hessian):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the BFGS-method

    Parameters
    ----------
    minimization_problem : minimization_problem
        An object containing information on the problem which we want to find
        a minimum
    xk : N-Array
        An array containing the current coordinates of our Newton method
    xk_1 : N-array
        An array containing the previous coordinates of our Newton method

    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the BFGS-method

    '''
    pass
        