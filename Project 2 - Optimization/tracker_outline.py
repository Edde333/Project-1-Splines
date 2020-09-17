# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:38:12 2020

@author: joelw
"""







"""
    Outline for possible tracker feature.
    You send in which parameters you want to have keep track of their 
    progression.
"""

class cc:
    
    
    def __init__(self, tracker):
        self.a = 12
        self.b = 4
        self.c = "aa"
        self.tracker = {t: [] for t in tracker}
    
    def fun(self):
        for i in range(10):
            a = i*5
            b = i**2
            c = str(i)
            for t in self.tracker:
                self.tracker[t].append(locals()[t])
        
        return self.a, self.tracker
            
        

if __name__ == '__main__':
    track = ['a','c']
    c = cc(track)
    
    out, tracked = c.fun()



