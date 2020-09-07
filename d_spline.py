# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
"""
This class is used to calculate a spline using a blossom algorithm. 
"""


class d_spline:
    
    
    def __init__(self, k, u, d):
        """
        Creates an object which is used to calculate the spline values.
        Parameters
        ----------
        k : Int
            The dimension of the spline (for example 1D, 2D, 3D...)
        u : Array (1xK-4)
            The chosen u-values of the spline
        d : Array (1xK-4)
            The chosen control points of the spline

        Returns
        -------
        None.

        """
        self.k = 3
        self.u = u
        self.d = d
       
    
    def dSpline(self,x):
        """
        The main function used to calculate the values of the spline
        
        Parameters
        ----------
        x : Float
            The points of interest of the spline. 

        Returns
        -------
        An array (1xn), or a matrix (kxn) depending on the 
            DESCRIPTION.
        """
        self.u = self.expandArray(self.u)
        ret = np.zeros((self.d.ndim,x.shape[0]))
        I = self.findHotInterval(self.u,x)
        # Run dSpline for each row(dimension) of the d matrix.
        if self.d.ndim > 1:
            for dimension in range(self.d.ndim):
                currd = self.expandArray(self.d[dimension])
                left = I
                right = I+1
                depth = 0
                ret[dimension] = self.dRecursive(x, I, currd, depth, left, right)
        else: 
            self.d = self.expandArray(self.d)
            left = I
            right = I+1
            depth = 0
            return self.dRecursive(x, I, self.d, depth, left, right)
        return ret


    def expandArray(self,u):
        """
        Expands the vector u by prepending k-1 elements
        equal to the first element in u, and appends 
        k-1 elements equal to the last element in u.
        
        Parameters
        ----------
        u : Array of floats (1,K-4)
            The vector that will be expanded

        Returns
        -------
        u : Array of floats (1,K)
            The extended Array u.

        """
        u = np.insert(u, 0, np.ones(self.k-1)*u[0])
        u = np.append(u,np.ones(self.k-1)*u[-1])
        return u
    
    
    def findHotInterval(self,u,x):
        """
        Calculates the hot interval for the point x.

        Parameters
        ----------
        u : Array of floats (1xK)
            The array of knots u_i.
        x : Array of floats (1xK)
            The array of points x_i to be evaluated

        Returns
        -------
        I : Array of ints.
            The array I contains the indices of the
            hot interval in the array u containing
            the knots.

        """
        I = np.zeros(len(x), dtype = 'int8')
        for i in range(len(x)):
            index = np.argmax(u >= x[i]) - 1
            if index < self.k-1: index = self.k-1
            I[i] = index
        return I
    
    
    def dRecursive(self, x, I, d, depth, left, right):
        """
        The recursive step of the blossom algorithm
        
        Parameters
        ----------
        x : Array of floats (1xK)
            The array of points x_i to be evaluated
        I : Array of ints.
            The array I contains the indices of the
            hot interval in the array u containing
            the knots.
        d : Array of floats
            An array of control points for the spline
        depth : Int
            Indicates at which depth of the recursion
            the current function is at
        left : Int
            The index of the left-most knot in the 
            interval I
        right : Int
            The index of the right-most knot in the
            interval I
        Returns
        -------
        The left-most value if in the deepest layer of recursion, else
        the calculated value for the given values in x.
        """
        if np.any(right-left == self.k+1):
            return d[left+1]
        a = self.alpha(x,left, right)
        d0 = self.dRecursive(x, I, d,  depth + 1, left - 1, right)
        d1 = self.dRecursive(x, I, d, depth + 1, left, right + 1)
        return a*d0 + (1-a)*d1
        
    
    def alpha(self, x,left, right):
        """
        Generates the values of alpha which are used to 
        calculate the values of x in the blossom algorithm.

        Parameters
        ----------
        x : Array of floats (1xK)
            The array of points x_i to be evaluated
       left : Int
            The index of the left-most knot in the 
            interval I
        right : Int
            The index of the right-most knot in the
            interval I

        Returns
        -------
        alph : FLoat
            Returns the value of alpha which is used to generate
            the values of x in the blossom algorithm

        """
        alphdivider = (self.u[right]-self.u[left])
        zeroValues = np.where(alphdivider == 0)
        alphdivider[zeroValues] = 1
        alph = (self.u[right] - x)/alphdivider
        alph[zeroValues] = 0
        return alph