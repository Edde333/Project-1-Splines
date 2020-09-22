# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:30:10 2020

@author: joelw
"""
import math

class cc:
    
    
    def __init__(self):
        self.a = 111
        
        
        
    def fun(self, tracking):
        track = {t : [] for t in tracking}
        x = 34.
        for i in range(10):
            x = x-2
            for name,l in track.items():
                l.append(locals()[name])
            
        return x, track
        


c = cc()
t = ['x']
x,tr = c.fun(t)

print()










