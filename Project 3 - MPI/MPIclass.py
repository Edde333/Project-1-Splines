#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 08:28:46 2020

@author: oliviaarnfalk
"""

import scipy as sp
from scipy import linalg
import numpy as np
from mpi4py import MPI 




class heatEqSolver():
    
    dx = 1/3
    omega = 0.8
    nbriterations = 10
    
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
        [temps, coords] = solver(self, p, dx, e_type, e_Temp, omega)
        comm.send(temps, 0)
        comm.send(coords, 0)
        
            
    
        
    def solver(self, p, dx, e_type, e_Temp, nbriterations, omega):
        
        xp = p[:,0]
        yp = p[:,1]
        
        x_max = np.argmax(xp)
        x_min = np.argmin(xp)
        y_max = np.argmax(yp)
        y_min = np.argmin(yp)
        
        width =  x_max - x_min
        height = y_max - y_min

        xx = np.linspace(0,width, width/dx)
        yy = np.linspace(0,height, height/dx)
        
        
        nodes = np.zeros((len(xx)*len(yy),3))
        a = 0
        for i in range (len(xx)) :
            for j in range (len(yy)):
                nodes[a,0] = xx[i]
                nodes[a,1] = yy[j]
                nodes[a,2] = None 
                a = a+1        
        
        coordinates = nodes + np.full((len(xx)*len(yy),3), [np.argmin(xp), np.argmin(yp), None])
        
        A = 1/dx**2 * get_diag_matrix(height/dx, len(v))
        
        
        x_dof = int(round(width / dx + 1))
        self.x_dof = x_dof
        y_dof = int(round(height / dx + 1))
        self.y_dof = y_dof
        nbr_dof = x_dof * y_dof
        self.nbr_dof = nbr_dof
        
        # Find which dof:s belong to the edges
        edof = [np.array([])] * len(p)
        for i in range(len(edof)):
            
            # Edge nbr. i-1
            e = np.array([p[i-1], p[i]])

            # Horisontal edge 
            if e[0,1] == e[1,1]:

                # Top edge (i = 0 in v_ij)
                if e[0,1] == y_max:
                    poss_dofs = range(x_dof)
                    # Loop through relevant dofs
                    for dof in poss_dofs:
                        if x[dof] >= e[0,0] and x[dof] < e[1,0]:
                            edof[i-1] = np.append(edof[i-1], dof)
                            
                # Bottom edge (i = y_dof in v_ij)
                else:
                    poss_dofs = range(x_dof*(y_dof-1), nbr_dof)
                    for dof in poss_dofs:
                        if x[dof] <= e[0,0] and x[dof] > e[1,0]:
                            edof[i-1] = np.append(edof[i-1], dof)

            # Vertical edge
            else:

                # Left edge (j = 0 in v_ij)
                if e[0,0] == x_min:
                    poss_dofs = range(0, nbr_dof, x_dof)
                    for dof in poss_dofs:
                        if y[dof] >= e[0,1] and y[dof] < e[1,1]:
                            edof[i-1] = np.append(edof[i-1], dof)

                # Right edge    
                else:
                    poss_dofs = range(x_dof-1, nbr_dof, x_dof)
                    for dof in poss_dofs:
                        if y[dof] <= e[0,1] and y[dof] > e[1,1]:
                            edof[i-1] = np.append(edof[i-1], dof)
                            
        self.edof = edof
        print(edof)

        
        
        
        
        # Get the A-matrix for the sub-problem
        A = 1/dx**2 * get_diag_matrix(height/dx, len(v))
        for i in range (len(nodes)):
            

        
        # Get the fbc vector of known boundaris with dirichlet conditions
        fbc = np.zeros(len(A)).T
        # Pairs boundaries to specific nodes, creates the constant vector fbc
        for i in range (len(e_type)):
            if e_type[i] == d:
                fbc = fbc + 1/dx**2 * T[i] @ indices
                # inserts:  - 1/dx**2 * T   into the f-vector
                   
                

        
        if rank != 1 :
            gammanodes = np.array([])
            gammaDir = 15*np.ones(len(gammanodes))
            
            
            for i in range (nbriterations):
                comm.send(gammaDir, src)
                gammaNeu = comm.recv(gammaNeu, source =1)
                # Update the fNbc vector with the values from rank == 1
                fNc = - 1/dx * gammaNeu #at the given instances 
                f = fbc + fNc
                vk1 = sp.linalg.solve(A,f)
                if i > 0:
                    vk1 = omega * vk1 + (1- omega) * vk
                vk = vk1
                # gamma takes the vk values calculated above that corresponds 
                # to the nodes along the boundary 
                gammaNeu = vk @ indices 
                
                
            temperatures = vk
            return [vk, coordinates]
                
                
                
            
        elif rank == 1:
            for n in range (nbriterations):
                f = fbc
                
                for i in range (len(src)):
                    gammaDir = comm.recv(gammaDir, source = src[i+2])
                    # add gamma to f-vector
                    fDc = 1/dx**2 * gammaDir @ indices
                    f = f + fDc
                
                vk = sp.linalg.solve(A,f)
                
                for i in range (len(src)):  
                    gammaNeu = vk @ indices 
                    comm.send(gammaNeu, src)
                    
            return [vk, coordinates]
                    
                    
                    
                    
            
        # Pairs the temperatures with the initial coordinates      
        # temperatures = 
            
        

if __name__ == '__main__'():
    points = np.array([[1,2], [2,2], [3,2], [3,1], [2,1], [2,0], [1,0], [0,0], 
                      [0,1], [1,1]]) 
    
    e_type_1 = [d,d,d,d,d,d]
    e_type_2 = [d,d,n,d]
    e_type_3 = [d,d,d,n]
    
    outer_edges = np.array([40, 40, 40, 15, 15, 5, 5, 40, 15, 15])
    
    print(points)
    
