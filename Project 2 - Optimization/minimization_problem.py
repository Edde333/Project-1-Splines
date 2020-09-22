import numpy as np

class minimization_problem():
    
    def __init__(self, function, guess, gradient = None):
        """

        Parameters
        ----------
        function : Python function that takes an N-array as input and returns a
            real number
        guess : N-array
            An array containing the coordinates of the inital guess of where the
            minimum is
        gradient : N-array of R^N to R functions, optional
            The gradient of the function for which we want to find a minimum.
            The default is None.

        Raises
        ------
        ValueError
            Raises a ValueError if the guess is not compatible with the function
            or the gradient

        Returns
        -------
        None.

        """
        
        self.function = function
        self.guess = guess
        self.dim = len(guess)
        
        # Checks that guess is applicable to the function (enough arguments)
        # Is it possible to check if there are too many arguments??
        try:
            function(guess)
        except ValueError:
            raise ValueError("'Guess' is not a compatible input to the function")
            
        # Creates gradient if not existant
        if gradient == None:
            self.gradient = lambda x: self.get_gradient(x)
        else:
            # Checks that guess is applicable to the gradient
            try:
                gradient(guess)
                self.gradient = gradient
            except ValueError:
                raise ValueError("'Guess' is not a compatible input to the gradient")  
            
    # For changing guess
    def __call__(self, guess):
        self.guess = guess
    
    def get_gradient(self, x):
        dx = 1e-6
        gradient = np.zeros(self.dim)
        
        # array used for difference calculation
        pc = np.identity(self.dim) * dx
        for i in range(self.dim):
            gradient[i] = (self.function(x+pc[i]) - self.function(x-pc[i])) / (2*dx)
        
        return gradient
    
if __name__ == "__main__":
    f = lambda x: x[0]**3 + x[1]**3 + x[2]**3
    guess = np.array([1,2,-3])
    mp = minimization_problem(f, guess)
    print(mp.function(guess))
    print(mp.gradient(guess))
