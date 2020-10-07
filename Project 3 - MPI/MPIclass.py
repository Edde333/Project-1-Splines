#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 08:28:46 2020

@author: oliviaarnfalk
"""

import scipy as sp
import numpy as np
from mpi4py import MPI 




class heatEqSolver():
    
    comm = MPI.COMM_WORLD()
    rank = comm.Get_rank()
    size = comm.Get_size()


    if rank == 0:
        for i in range(size-1):
            temps = comm.recv(source = i + 1)
            coords = comm.recv(source = i + 1)
    # plots the  temperature distribution       
            
            
            
    else:
        # find sub-problem, use solver to get temperatures
        [temps, coords] = solver()
        comm.send(temps, 0)
        comm.send(coords, 0)
        
            
    
        
    def solver(p, dx, e_type, e_Temp, nbriterations):
        
        xp = p[:,0]
        yp = p[:,1]
        
        width = np.argmax(xp) - np.argmin(xp)
        height = np. argmax(yp) - np.argmin(yp)

        xx = np.linspace(0,width, width/dx)
        yy = np.linspace(0,height, height/dx)
        
        
        nodes = np.zeros((len(xx)*len(yy),2))
        
        
        a = 0
        for i in range (len(xx)) :
            for j in range (len(yy)):
                nodes[a,0] = xx[i]
                nodes[a,1] = yy[j]
                a = a+1
        
        #coordinates = 
        
        
        # Get the A-matrix for the sub-problem
        A = get_diag_matrix(height/dx, len(v))
        
        
        # Pairs boundaries to specific nodes, creates the constant vector f
        for i in range (len(e_type)):
            if e_type[i] == d:
                # inserts:  - 1/dx**2 * T   into the f-vector, 
                 
        
        
        # Get the f vector of known temperatures with dirichlet conditions
        f = np.zeros(len(nodes)).T
        
        
        
        
        
        fetch = np.array([])
        e_initint = np.array([])
        edge_nodes = np.array([])
        src = np.array([])

        
        if rank != 1 :
            gammanodes = np.array([])
            gamma{rank} = 15*np.ones(len(gammanodes))
            
            comm.send(edgetemp, src)
            
            for i in range (nbriterations):
                gammaNeu = comm.recv(gammaNeu, source =1)
                
                
            
        elif rank == 1:
            for n in range (nbriterations)
                for i in range (len(src))
                    gamma{src[i]} = comm.recv(gamma{i}, source = src{i})
                    # add gamma to f-vector
                    
                    v = sp.linalg.solve(A,f)
                    
            
        # Pairs the temperatures with the initial coordinates      
        temps =     
            
        

if __name__ == '__main__'():
    points = np.array([[1,2], [2,2], [3,2], [3,1], [2,1], [2,0], [1,0], [0,0], 
                      [0,1], [1,1]]) 
    
    e_type_1 = [d,d,d,d,d,d]
    e_type_2 = [d,d,n,d]
    e_type_3 = [d,d,d,n]
    
    outer_edges = np.array([40, 40, 40, 15, 15, 5, 5, 40, 15, 15])
    
    print(points)
    
