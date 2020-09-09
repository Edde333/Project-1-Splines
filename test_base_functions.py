
import unittest
from  scipy import *
from  pylab import *
from base_functions import *
from d_spline import *
import matplotlib.pyplot as plt
import numpy as np

class TestIdentity(unittest.TestCase): 
    

    
    def test_baseFuncSumEqualOne (self):
        x = linspace(0,1,10)
        u = linspace(0,1,26)
        u[ 1] = u[ 2] = u[ 0]
        u[-3] = u[-2] = u[-1]
        result = baseFuncSum(x, u)
        expected = 1
        for i in range(10):
            self.assertAlmostEqual(result[i] , expected)



    def test_BaseFunc (self):
        x = [0.2]
        u = linspace(0,1,26)
        u[ 1] = u[ 2] = u[ 0]
        u[-3] = u[-2] = u[-1]
        d = [(-12.73564, 9.03455),
         (-26.77725, 15.89208),
         (-42.12487, 20.57261),
         (-15.34799, 4.57169),
         (-31.72987, 6.85753),
         (-49.14568, 6.85754),
         (-38.09753, -1e-05),
         (-67.92234, -11.10268),
         (-89.47453, -33.30804),
         (-21.44344, -22.31416),
         (-32.16513, -53.33632),
         (-32.16511, -93.06657),
         (-2e-05, -39.83887),
         (10.72167, -70.86103),
         (32.16511, -93.06658),
         (21.55219, -22.31397),
         (51.377, -33.47106),
         (89.47453, -33.47131),
         (15.89191, 0.00025),
         (30.9676, 1.95954),
         (45.22709, 5.87789),
         (14.36797, 3.91883),
         (27.59321, 9.68786),
         (39.67575, 17.30712)]
        result = getCubicSpline(x,u,d)
        a = result[0]
        expectedz=  -31.90219167
        expectedy = 6.47655833
        self.assertAlmostEqual(a[0] , expectedz)
        self.assertAlmostEqual(a[1], expectedy)
        


    def test_xOutsideInterval(self):
        u = linspace(0,1,26)
        u[ 1] = u[ 2] = u[ 0]
        u[-3] = u[-2] = u[-1]
        d = [(-12.73564, 9.03455),
         (-26.77725, 15.89208),
         (-42.12487, 20.57261),
         (-15.34799, 4.57169),
         (-31.72987, 6.85753),
         (-49.14568, 6.85754),
         (-38.09753, -1e-05),
         (-67.92234, -11.10268),
         (-89.47453, -33.30804),
         (-21.44344, -22.31416),
         (-32.16513, -53.33632),
         (-32.16511, -93.06657),
         (-2e-05, -39.83887),
         (10.72167, -70.86103),
         (32.16511, -93.06658),
         (21.55219, -22.31397),
         (51.377, -33.47106),
         (89.47453, -33.47131),
         (15.89191, 0.00025),
         (30.9676, 1.95954),
         (45.22709, 5.87789),
         (14.36797, 3.91883),
         (27.59321, 9.68786),
         (39.67575, 17.30712)]
        x = [ 2]
        
        with self.assertRaises(Exception): getCubicSpline(x,u,d)

    
    def testL2Norm(k = 3, u = np.linspace(1,10,10), x = linspace(0,10,100),
               d = np.array([[0, 1, 3, 3, 4, 5, 6, 7, 8, 9, 10],
              [2, 4, 2, 7, 4, 8, 5, 3, 4, 6, 0]])):
        k = 3
        u = np.linspace(1,10,10)
        u = np.insert(u, 0, np.ones(k-1)*u[0])
        u = np.append(u2,np.ones(k-1)*u[-1])
        x = linspace(1,10,100)
        d = np.array([[1, 3, 3, 4, 5, 6, 7, 8, 9, 10],
                  [4, 2, 7, 4, 8, 5, 3, 4, 6, 0]])
        d1 = np.array([[1,1,1, 3, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10],
                  [4, 4, 4, 2, 7, 4, 8, 5, 3, 4, 6, 0,0,0]])
    
        ds = d_spline(k, u, d)
        res = ds.dSpline(x,3)
        

        d2x = np.insert(d[0], 0, np.ones(k-1)*d[0][0])
        d2y = np.insert(d[1], 0, np.ones(k-1)*d[1][0])
        N_i = []
        for i in range(len(u2)-2):
            N_i.append(getBaseFunc(u2,i))
            
        base_x = N_i[0](x)*d2x[0]
        base_y = N_i[0](x)*d2y[0]
    
        for i in range(len(N_i)-1):
            base_x +=  N_i[i+1](x)*d2x[i+1]
            base_y +=  N_i[i+1](x)*d2y[i+1]
            
        L2norm = max(sqrt((res[0]-base_x)**2 + (res[1] - base_y)**2))
        expected = 0
        
        self.assertAlmostEqual(L2norm, expected)
        
        """
        plt.figure()
        plt.plot(res[0],res[1])
        plt.show()
        plt.scatter(d[0],d[1])
        plt.figure()
        plt.plot(base_x,base_y)
        plt.show()
        plt.scatter(d2x,d2y)
        """   
        
    
if __name__ == "__main__":
    testL2Norm()

