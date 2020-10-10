import numpy as np
from get_mesh import get_mesh
from scipy.sparse import dia_matrix
from scipy import linalg

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
      
        # Initiate A
        ex = np.ones(self.nbr_dof)  
        data = np.array([ex,-4*ex,ex])
        offsets = np.array([-1, 0,  1])
        A = dia_matrix((data, offsets), shape=(self.nbr_dof, self.nbr_dof)).toarray()      
        sub1 = np.eye(self.nbr_dof, k = -self.x_dof)
        sub2 = np.eye(self.nbr_dof, k = self.x_dof)   
        A = (A + sub1 + sub2)
        
        # Initiate f
        f = np.zeros(self.nbr_dof) 
        
        # Modify A and f by looping through edge_dofs
        for i in range(len(self.edof)):
            # Edge nbr. i
            
            # Edge orientation
            eo = self.edge_orient[i]
            
            # Loop through dofs in edge nbr. i
            for dof in self.edof[i]:
                
                # Diagonal element 
                A[dof,dof] = -3
                
                # Dirichlet condition (No A modification because corresponding
                # row and column will be reduced away)
                if self.edge_type[i] == 'd':
                    # Condition is known and constant (constructs constant f)
                    if self.fetch[i] == None:
                        fv = (-1) * self.edge_init[i] / (self.dx**2)
                        if eo == 't':
                            f[dof+self.x_dof] += fv
                            if dof != 0:    
                                f[dof-1] += fv
                            if dof != self.x_dof-1:
                                f[dof+1] += fv
                        elif eo == 'b':
                            f[dof-self.x_dof] += fv
                            if dof != self.nbr_dof-self.x_dof:
                                f[dof-1] += fv
                            if dof != self.nbr_dof-1:
                                f[dof+1] += fv
                        elif eo == 'l':
                            f[dof+1] += fv
                            if dof != 0:
                                f[dof-self.x_dof] += fv
                            if dof != self.nbr_dof-self.x_dof:
                                f[dof+self.x_dof] += fv
                        else: # eo == 'r'
                            f[dof-1] += fv
                            if dof != self.x_dof-1:
                                f[dof-self.x_dof] += fv
                            if dof != self.nbr_dof-1:
                                f[dof+self.x_dof] += fv
                        
                
                # Neumann condition (A is modified and basic f constructed)
                else:
                    # A modification
                    if eo == 'l' and dof != 0:
                        A[dof, dof-1] = 0
                    elif eo == 'r' and dof != self.nbr_dof-1:
                        A[dof, dof+1] = 0
                    # eo == 't' and 'b' is taken care of automatically
                    
                    # Condition is known and constant
                    if self.fetch[i] == None:  
                        f[dof] += (-1) * self.edge_init[i] / self.dx
        
        # Handle vertices seperately        
        v_dof = np.array([0, self.x_dof-1, 
                          self.nbr_dof-self.x_dof, self.nbr_dof-1])
        # If Neumann node
        if not (v_dof[0] in self.d_dofs) and (v_dof[0]+1 in self.d_dofs):
            A[v_dof[0], v_dof[0]+1] = 0       
        if (not v_dof[1] in self.d_dofs) and (v_dof[1]-1 in self.d_dofs):
            A[v_dof[1], v_dof[1]-1] = 0   
        if (not v_dof[2] in self.d_dofs) and (v_dof[2]+1 in self.d_dofs):
            A[v_dof[2], v_dof[2]+1] = 0
        if (not v_dof[3] in self.d_dofs) and (v_dof[3]-1 in self.d_dofs):
            A[v_dof[3], v_dof[3]-1] = 0
                
                
        # Complete A definition  
        A = A / (self.dx**2)
        # Define reduced A-matrix
        A_red = A[self.red_dofs][:,self.red_dofs]
        
        # Add known bc:s to v-vector
        v_old = np.zeros(self.nbr_dof)
        v = np.zeros(self.nbr_dof)
        for i in range(len(self.edof)):
            if self.fetch[i] == None and self.edge_type[i] == 'd':
                v[self.edof[i]] = self.edge_init[i]
        
        # Solving loop        
        for i in range(nbr_iter):
            
            # Recieve fetched conditions and update f
            f_new = f.copy()
            # Loop through edges
            for j in range(len(self.edof)):
                # Edge nbr. j
                
                # Edge that fetches boundary condition from other process
                if self.fetch[j] != None:
                    cond = comm.recv(source = self.fetch[j])

                    eo = self.edge_orient[j]
                    
                    # Dirichlet condition
                    if self.edge_type[j] == 'd':
                        v[self.edof[j]] = cond
                        
                        fv = (-1) * cond / (self.dx**2)
                        
                        if eo == 't':
                            f_new[self.edof[j]+self.x_dof] += fv
                        elif eo == 'b':
                            f_new[self.edof[j]-self.x_dof] += fv
                        elif eo == 'l':
                            f_new[self.edof[j]+1] += fv
                        else: # eo == 'r'
                            f_new[self.edof[j]-1] += fv
                           
                    # Neumann condition
                    else:
                        f_new[self.edof[j]] += (-1) * cond / self.dx
            
            
            # Reduce equation system and solve
            f_red = f_new[self.red_dofs]
            v_red = linalg.inv(A_red).dot(f_red)
            v[self.red_dofs] = v_red

            # Relaxation
            if i > 0:
                v = omega * v + (1-omega) * v_old
            v_old = v.copy()
            
            # Loop through edges
            for j in range(len(self.edof)):
                # Edge nbr. j
                
                # Edge that sends boundary condition from other process
                if self.fetch[j] != None:
                    # Send Dirichlet
                    if self.edge_type[j] == 'n':
                        data = v[self.edof[j]]
                    # Send Neumann    
                    else:
                        eo = self.edge_orient[j]
                        if eo == 'l':
                            data = (v[self.edof[j]+1] - v[self.edof[j]]) / self.dx
                        elif eo == 'r':
                            data = (v[self.edof[j]-1] - v[self.edof[j]]) / self.dx
                        elif eo == 't':
                            data = (v[self.edof[j]+self.x_dof] - v[self.edof[j]]) / self.dx
                        else: # eo == 'b'
                            data = (v[self.edof[j]-self.x_dof] - v[self.edof[j]]) / self.dx

                    # Don't send last time if in process 2 or 3
                    if i == nbr_iter-1:
                        if comm.Get_rank() == 1:
                            comm.send(data, self.fetch[j])
                    else:    
                        comm.send(data, self.fetch[j])
   
        return np.reshape(v, (int(self.y_dof),int(self.x_dof)))

if __name__ == "__main__":
    guess = 20
    dx = 0.15
    # Define region
    points = np.array([(0,2), (1,2), (1,1), (1,0), (0,0), (0,1)])
    #points = np.array([(0,1), (1,1), (1,0), (0,0)])
    edge_type = np.array(['d', 'd', 'd', 'd', 'd', 'd'])
    fetch = np.array([None, 3, None, None, 2, None])
    edge_init = np.array([40, None, 5, 15, None, 40])

    # Create region
    r = region(points, edge_type, fetch, edge_init, dx)
    
    # Solve
    #v = r.solve(1)