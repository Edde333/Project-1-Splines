# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 10:42:22 2020

@author: joelw
"""



def fun(x,y,z, a=23,b=4,c=5):
    
    return x+y+z+a+b+c
    
    


f = lambda x: fun(*x)

q = 12
w = -111
e = -4.5
input_log = ["q", "w", "e"]
input_var = [globals()[i] for i in input_log]
print( f((10,23,44,33)) )







