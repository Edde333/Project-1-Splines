# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 10:42:22 2020

@author: joelw
"""



def fun(x,y,z, a=None,b=None,c=None):
    input_var = locals()
    for name,value in input_var.items():
        if value is None: 
            input_var[name] = 123.
    print(locals())
    locals()['a'] = "ASSADDDSASADSDSDSDA"
    print(locals())
    return x+y+z+a+b+c
    
    


f = lambda x: fun(*x)

q = 12
w = -111
e = -4.5
input_log = ["q", "w", "e"]
input_var = [globals()[i] for i in input_log]
fun(1,2,2,b=3)








