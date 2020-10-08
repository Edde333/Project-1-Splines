#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:16:37 2020

@author: oliviaarnfalk
"""

import scipy as sp
from scipy import linalg
import numpy as np
from mpi4py import MPI 
from Diag_matrix import *

def Create_matrix(self, p, edge_type, fetch, edge_init, dx):
        xp = p[:,0]
        yp = p[:,1]

        x_max = max(xp)
        x_min = min(xp)
        y_max = max(yp)
        y_min = min(yp)
        

        width =  x_max - x_min
        height = y_max - y_min
        
        xx = np.linspace(0, width, int(width/dx))
        yy = np.linspace(0, height, int(height/dx))
        
        # Create vector for all teh nodes
        nodes = np.zeros((len(xx)*len(yy),2))
        a = 0
        for i in range (len(xx)) :
            for j in range (len(yy)):
                nodes[a,0] = xx[i]
                nodes[a,1] = yy[j]
                a = a+1  
                
        print(nodes)
        coordinates = nodes + np.full((len(xx)*len(yy),2), [xmin, xmax])
        print(coordinates)
        
        A = 1/dx**2 * get_diag_matrix(height/dx, len(nodes))
        
        
        
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
        
        
        print(edof[1])
        edof_Dirichlet = np.array([])
        edof_Neumann = np.array([])
        for i in range (e_type):
            if e_type[i] == d:
                edof_Dirichlet = np.append(edof_Dirichlet, edof[i])
            else: 
                edof_Neumann = np.append(edof_Neumann, edof[i])
          
            
        
        Areduced = A [[notdiricehlet],[notdirichlet]]
        
        
        
        
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
              
                
if __name__ == "__main__":
    guess = 20
    dx = 0.25
    # Define region
    points = np.array([(0,2), (1,2), (1,1), (1,0), (0,0), (0,1)])
    
    edge_type = np.array(['d', 'd', 'd', 'd', 'd', 'd'])
    fetch = np.array([None, 3, None, None, None, None])
    edge_init = np.array([40, guess, 15, 5, guess, 15])

    # Create region
    r = Create_matrix(points, edge_type, fetch, edge_init, dx)
    