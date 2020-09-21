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
    elif hessian_approximation_method == "finite_differences":
        return finite_differences(gradient, xk)
    elif hessian_approximation_method == "brute_inverse_hessian":
        return brute_inv_hessian(function, xk)
    else:
        output = "There is no method called: " + hessian_approximation_method
        raise ValueError(output)



def finite_differences(gradient, x):
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
    deltak = xk - xk_1
    gammak = gradient(xk) - gradient(xk_1)
    u = deltak - prev_hessian@gammak
    a = 1/(u@gammak)
    new_hess = prev_hessian + a*np.outer(u,u)
    return new_hess

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
    deltak = xk - xk_1
    gammak = gradient(xk) - gradient(xk_1)

    first_term = (np.outer(deltak,deltak))/(np.inner(deltak,gammak))
    sec_term = (prev_hessian@np.outer(gammak,gammak@prev_hessian))/(gammak@prev_hessian@gammak)

    return prev_hessian + first_term - sec_term

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
    deltak = xk - xk_1
    gammak = gradient(xk) - gradient(xk_1)
    
    first_term = 1 + (np.inner(gammak,prev_hessian@gammak))/(np.inner(deltak,gammak))
    sec_term = (np.outer(deltak,deltak))/(np.inner(deltak,gammak))
    third_term = (np.outer(deltak,gammak@prev_hessian) + prev_hessian@np.outer(gammak,deltak))/(np.inner(deltak,gammak))

    return prev_hessian + first_term*sec_term - third_term 


def brute_inv_hessian(function, xk, h = 0.001):
    '''
    Calculates the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, using a simple brute force method
    Parameters
    ----------
    function : Lamda function
        The function that should be minimized.
    xk : N-Array
        An array containing the current coordinates of our Newton method
    Returns
    -------
    Returns the approximated inverse Hessian matrix for the probmlem at the
        specified points xk, and xk_1 using the BFGS-method
    '''
    input_dim = xk.size
    hessian = np.empty( ( input_dim, input_dim))
    I =np.zeros(input_dim)
    J = np.zeros( input_dim)
    for index in np.ndindex(hessian.shape):
        I[index[0]] = h
        J[index[1]] = h
        hessian[index] = (function(xk+I+J) + function(xk-I-J) - function(xk-I+J) - function(xk+I-J))/(4*h*h)
        I[index[0]] = 0
        J[index[1]] = 0
    hessian = 1/2*(hessian + hessian.T)
        # checks that hessian is positive definite
    try:
        np.linalg.cholesky(hessian)
    except np.linalg.LinAlgError:
        raise Exception("Hessian matrix not positive definite")
    inv_hessian = np.linalg.inv(hessian)
    return inv_hessian



if __name__ == "__main__":
    function = lambda x: x[0]**3 + x[1]**3 + x[2]**3
    gradient = lambda x: np.array([6*x[0],6*x[1],6*x[2]])
    xk_1 = np.array([1, 1, 1])
    xk = np.array([2, 2, 2])
    prev_hessian = np.array([[6, 0, 0],[0,6, 0],[0,0,6]])
    y = function(xk)
    new_hessian = BFGS(function,gradient, xk, xk_1, prev_hessian)
    print(new_hessian)
