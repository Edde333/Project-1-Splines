import numpy as np
from get_mesh import get_mesh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

class region:
    """
    Rectangular region.

    """

    def __init__(self, points, edge_type, fetch, edge_init, dx):
        
        # Save attributes
        self.edge_type = edge_type
        self.fetch = fetch
        self.edge_init = edge_init
        self.dx = dx
        
        # Create mesh
        x_dof, y_dof, edof, edge_orient = get_mesh(points, dx, fetch)
        nbr_dof = x_dof * y_dof
        
        self.x_dof = x_dof
        self.y_dof = y_dof
        self.nbr_dof = nbr_dof
        self.edof = edof
        self.edge_orient = edge_orient
        
        # Find the dirichlet dofs
        d_dofs = np.array([], dtype = int)
        for i in range(len(edof)):
            if edge_type[i] == 'd':
                d_dofs = np.append(d_dofs, edof[i])
        
        # Find the dofs to keep in reduced equation system
        red_dofs = np.array([], dtype = int)
        for i in range(nbr_dof):
            if not (i in d_dofs):
                red_dofs = np.append(red_dofs, i)
        
        self.d_dofs = d_dofs
        self.red_dofs = red_dofs
        

    def solve(self, comm, nbr_iter = 10, omega = 0.8):
        """
        Solves the Laplace equation for the region and sends data to
        corresponding process through MPI.

        Returns
        -------
        Temperature field, (x_dof * y_dof)-np-array

        """
        
        # Initiate sparse matrix representation of equation system
        values = np.array([])
        row_ind = np.array([], dtype = int)
        col_ind = np.array([], dtype = int)
        
        f_values = np.array([])
        f_ind = np.array([], dtype = int)
        
        # A-matrix
        # Diagonal
        values = np.append(values, np.ones(self.nbr_dof) * (-4))
        row_ind = np.append(row_ind, np.arange(self.nbr_dof))
        col_ind = np.append(col_ind, np.arange(self.nbr_dof))
        
        # Diagonal offset 1
        values = np.append(values, np.ones(2 * self.nbr_dof - 2))
        # Top diagonal
        row_ind = np.append(row_ind, np.arange(self.nbr_dof-1))
        col_ind = np.append(col_ind, np.arange(1, self.nbr_dof))
        # Bottom diagonal
        row_ind = np.append(row_ind, np.arange(1, self.nbr_dof))
        col_ind = np.append(col_ind, np.arange(self.nbr_dof-1))
        
        # Diagonal offset x_dof
        values = np.append(values, np.ones(2*(self.nbr_dof-self.x_dof)))
        # Top diagonal
        row_ind = np.append(row_ind, np.arange(self.nbr_dof-self.x_dof))
        col_ind = np.append(col_ind, np.arange(self.x_dof, self.nbr_dof))
        # Bottom diagonal         
        row_ind = np.append(row_ind, np.arange(self.x_dof, self.nbr_dof))
        col_ind = np.append(col_ind, np.arange(self.nbr_dof-self.x_dof))
        
        
        # Modify A and f by looping through edge_dofs
        for i in range(len(self.edof)):
            # Edge nbr. i
            
            # Edge orientation
            eo = self.edge_orient[i]
            
            # Loop through dofs in edge nbr. i
            for dof in self.edof[i]:
                
                # Dirichlet condition (No A modification because corresponding
                # row and column will be reduced away)
                if self.edge_type[i] == 'd':
                    # Condition is known and constant (constructs constant f)
                    if self.fetch[i] == None:
                        f_values = np.append(f_values, (-1) * self.edge_init[i] / (self.dx**2))
                        
                        if eo == 't':
                            f_ind = np.append(f_ind, int(dof+self.x_dof))
                        elif eo == 'b':
                            f_ind = np.append(f_ind, int(dof-self.x_dof))
                        elif eo == 'l':
                            f_ind = np.append(f_ind, int(dof+1))
                        else: # eo == 'r'
                            f_ind = np.append(f_ind, int(dof-1))
                        
                
                # Neumann condition (A is modified and basic f constructed)
                else:
                    # A modification
                    values[dof] = -3
                    
                    if eo == 'l' and dof != 0:
                        values = np.append(values, -1)
                        row_ind = np.append(row_ind, dof)
                        col_ind = np.append(col_ind, dof-1)
                    elif eo == 'r' and dof != self.nbr_dof-1:
                        values = np.append(values, -1)
                        row_ind = np.append(row_ind, dof)
                        col_ind = np.append(col_ind, dof+1)
                    # eo == 't' and 'b' is taken care of automatically
                    
                    # Condition is known and constant
                    if self.fetch[i] == None:  
                        f_values = np.append(f_values, (-1) * self.edge_init[i] / self.dx)
                        f_ind = np.append(f_ind, dof)
        
        # Handle vertices seperately        
        v_dof = np.array([0, self.x_dof-1, 
                          self.nbr_dof-self.x_dof, self.nbr_dof-1],
                         dtype = int)
            
        for i in range(len(self.edof)):
            # Neumann edge
            if self.edge_type[i] == 'n':
                # Loop through the vertices
                for dof in v_dof:
                    if dof in self.edof[i]:
                        # Center diagonal
                        values = np.append(values, -1)
                        row_ind = np.append(row_ind, dof)
                        col_ind = np.append(col_ind, dof)
                        
                        # Top edge
                        if self.edge_type[i] == 't':
                            # Top right corner
                            if dof == v_dof[1]:
                                values = np.append(values, -1)
                                row_ind = np.append(row_ind, dof)
                                col_ind = np.append(col_ind, dof+1)
                                
                                adj_edge = i+1
                                
                            # Top left corner
                            else:
                                adj_edge = -1
                        
                        # Bottom edge
                        elif self.edge_type[i] == 'b':
                            # Bottom left corner
                            if dof == v_dof[2]:
                                values = np.append(values, -1)
                                row_ind = np.append(row_ind, dof)
                                col_ind = np.append(col_ind, dof-1)
                                
                                adj_edge = i+1
                            
                            # Bottom right corner
                            else:
                                adj_edge = i-1
                                    
                        # Left edge
                        elif self.edge_type[i] == 'l':
                            # Top left corner
                            if dof == v_dof[0]:    
                                adj_edge = 0
                            # Bottom left corner
                            else:
                                adj_edge = i-1
                                
                        # Right edge
                        else:
                            # Top right corner
                            if dof == v_dof[1]:
                                adj_edge = i-1
                            # Bottom right corner
                            else:
                                adj_edge = i+1
                                
                        # Adjacent dirichlet edge gives f value
                        if self.edge_type[adj_edge] == 'd':
                            # Other edge is not fetched since it would
                            # include the corner
                            f_values = np.append(f_values, (-1) * edge_init[adj_edge] / (self.dx**2))
                            f_ind = np.append(f_ind, dof)
                                
                
        # Complete A definition  
        values = values / (self.dx**2)
        
        # Define reduced A-matrix
        # Reduce rows
        red_row_ind = np.array([], dtype = int)
        red_col_ind = np.array([], dtype = int)
        red_values = np.array([])
        for i in range(len(self.red_dofs)):
            indices = np.where(row_ind == self.red_dofs[i])
            red_values = np.append(red_values, values[indices])
            red_row_ind = np.append(red_row_ind, np.ones(len(indices[0]), dtype = int) * i)
            red_col_ind = np.append(red_col_ind, col_ind[indices])
            
        values = red_values.copy()
        row_ind = red_row_ind.copy()
        col_ind = red_col_ind.copy()
        
        # Reduce columns
        red_row_ind = np.array([], dtype = int)
        red_col_ind = np.array([], dtype = int)
        red_values = np.array([])
        for i in range(len(self.red_dofs)):
            indices = np.where(col_ind == self.red_dofs[i])
            red_values = np.append(red_values, values[indices])
            red_row_ind = np.append(red_row_ind, row_ind[indices])
            red_col_ind = np.append(red_col_ind, np.ones(len(indices[0]), dtype = int) * i)
        
        A = csr_matrix((red_values, (red_row_ind, red_col_ind)))
        
        # # Add known bc:s to v-vector
        # v_old = np.zeros(self.nbr_dof)
        # v = np.zeros(self.nbr_dof)
        # for i in range(len(self.edof)):
        #     if self.fetch[i] == None and self.edge_type[i] == 'd':
        #         v[self.edof[i]] = self.edge_init[i]
        
        # # Solving loop        
        # for i in range(nbr_iter):
            
        #     # Recieve fetched conditions and update f
        #     f_new = f.copy()
        #     # Loop through edges
        #     for j in range(len(self.edof)):
        #         # Edge nbr. j
                
        #         # Edge that fetches boundary condition from other process
        #         if self.fetch[j] != None:
        #             cond = comm.recv(source = self.fetch[j])

        #             eo = self.edge_orient[j]
                    
        #             # Dirichlet condition
        #             if self.edge_type[j] == 'd':
        #                 v[self.edof[j]] = cond
                        
        #                 fv = (-1) * cond / (self.dx**2)
                        
        #                 if eo == 't':
        #                     f_new[self.edof[j]+self.x_dof] += fv
        #                 elif eo == 'b':
        #                     f_new[self.edof[j]-self.x_dof] += fv
        #                 elif eo == 'l':
        #                     f_new[self.edof[j]+1] += fv
        #                 else: # eo == 'r'
        #                     f_new[self.edof[j]-1] += fv
                           
        #             # Neumann condition
        #             else:
        #                 f_new[self.edof[j]] += (-1) * cond / self.dx
            
            
        #     # Reduce equation system and solve
        #     f_red = f_new[self.red_dofs]
        #     v_red = linalg.inv(A_red).dot(f_red)
        #     v[self.red_dofs] = v_red

        #     # Relaxation
        #     if i > 0:
        #         v = omega * v + (1-omega) * v_old
        #     v_old = v.copy()
            
        #     # Loop through edges
        #     for j in range(len(self.edof)):
        #         # Edge nbr. j
                
        #         # Edge that sends boundary condition from other process
        #         if self.fetch[j] != None:
        #             # Send Dirichlet
        #             if self.edge_type[j] == 'n':
        #                 data = v[self.edof[j]]
        #             # Send Neumann    
        #             else:
        #                 eo = self.edge_orient[j]
        #                 if eo == 'l':
        #                     data = (v[self.edof[j]+1] - v[self.edof[j]]) / self.dx
        #                 elif eo == 'r':
        #                     data = (v[self.edof[j]-1] - v[self.edof[j]]) / self.dx
        #                 elif eo == 't':
        #                     data = (v[self.edof[j]+self.x_dof] - v[self.edof[j]]) / self.dx
        #                 else: # eo == 'b'
        #                     data = (v[self.edof[j]-self.x_dof] - v[self.edof[j]]) / self.dx

        #             # Don't send last time if in process 2 or 3
        #             if i == nbr_iter-1:
        #                 if comm.Get_rank() == 1:
        #                     comm.send(data, self.fetch[j])
        #             else:    
        #                 comm.send(data, self.fetch[j])
   
        # return np.reshape(v, (int(self.y_dof),int(self.x_dof)))

if __name__ == "__main__":
    guess = 20
    dx = 0.01
    # Define region
    points = np.array([(0,2), (1,2), (1,1), (1,0), (0,0), (0,1)])
    #points = np.array([(0,1), (1,1), (1,0), (0,0)])
    edge_type = np.array(['d', 'd', 'd', 'd', 'd', 'd'])
    fetch = np.array([None, 3, None, None, 2, None])
    edge_init = np.array([40, None, 5, 15, None, 40])

    # Create region
    r = region(points, edge_type, fetch, edge_init, dx)
    
    # Solve
    r.solve(1)