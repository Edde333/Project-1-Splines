from mpi4py import MPI
import numpy as np
import sys
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
        x_min = sys.float_info.max
        y_min = sys.float_info.max

        x_max = (-1) * sys.float_info.max
        y_max = (-1) * sys.float_info.max
        
        for p in points:
            if p[0] < x_min:
                x_min = p[0]
            if p[0] > x_max:
                x_max = p[0]
            
            if p[1] < y_min:
                y_min = p[1]
            if p[1] > y_max:
                y_max = p[1]
        
        base = x_max - x_min
        height = y_max - y_min

        x_dof = int(round(base / dx + 1))
        self.x_dof = x_dof
        y_dof = int(round(height / dx + 1))
        self.y_dof = y_dof
        nbr_dof = x_dof * y_dof
        self.nbr_dof = nbr_dof
        
        # Coordinate vectors
        x = np.array([])
        xr = np.arange(x_dof) * dx
        for i in range(y_dof):
            x = np.append(x, xr)

        y = np.array([])
        yr = np.ones(x_dof) * height
        for i in range(y_dof):
            y = np.append(y, yr)
            yr -= dx

        # Find which dof:s belong to the edges
        edof = [np.array([], dtype = int)] * len(points)
        edge_orient = [None] * len(points)
        d_dofs = np.array([], dtype = int) # Dirichlet dofs
        for i in range(len(edof)):
            # Edge nbr. i-1
            e = np.array([points[i-1], points[i]])

            # Horisontal edge 
            if e[0,1] == e[1,1]:

                # Top edge (i = 0 in v_ij)
                if e[0,1] == y_max:
                    edge_orient[i-1] = 't'

                    poss_dofs = range(x_dof)
                    # Loop through relevant dofs
                    for dof in poss_dofs:
                        if x[dof] > e[0,0] and x[dof] < e[1,0]:
                            edof[i-1] = np.append(edof[i-1], dof)
                # Bottom edge (i = y_dof in v_ij)
                else:
                    edge_orient[i-1] = 'b'
                    
                    poss_dofs = range(x_dof*(y_dof-1), nbr_dof)
                    for dof in poss_dofs:
                        if x[dof] < e[0,0] and x[dof] > e[1,0]:
                            edof[i-1] = np.append(edof[i-1], dof)

            # Vertical edge
            else:

                # Left edge (j = 0 in v_ij)
                if e[0,0] == x_min:
                    edge_orient[i-1] = 'l'
                    
                    poss_dofs = range(0, nbr_dof, x_dof)
                    for dof in poss_dofs:
                        if y[dof] >= e[0,1] and y[dof] <= e[1,1]:
                            edof[i-1] = np.append(edof[i-1], dof)

                # Right edge    
                else:
                    edge_orient[i-1] = 'r'
                    
                    poss_dofs = range(x_dof-1, nbr_dof, x_dof)
                    for dof in poss_dofs:
                        if y[dof] <= e[0,1] and y[dof] >= e[1,1]:
                            edof[i-1] = np.append(edof[i-1], dof)
            
            # Keep the dofs sorted for correct communication
            edof[i-1].sort()

            # Save the dirichlet dofs
            if edge_type[i-1] == 'd':
                d_dofs = np.append(d_dofs, edof[i-1])
                
        # Handle mulpiple edges on vertical edge
        for i in range(len(edof)):
            if edge_orient[i-1] == edge_orient[i]:
                # If we are fetching we have dof priority
                if fetch[i-1] != None:
                    if edge_orient[i] == 'l':
                        edof[i] = np.delete(edof[i], -1)
                    else: # right edge
                        edof[i] = np.delete(edof[i], 0)
            
        # Find the dofs to keep in reduced equation system
        red_dofs = np.array([], dtype = int)
        for i in range(nbr_dof):
            if not (i in d_dofs):
                red_dofs = np.append(red_dofs, i)
                
                
        self.edof = edof
        self.edge_orient = edge_orient
        self.red_dofs = red_dofs

    def solve(self, comm):
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
        f = np.zeros(self.nbr_dof) 
        
        # Modify A and create basic f by looping through edge_dofs
        for i in range(len(self.edof)):
            # Edge nbr. i
            
            # Edge orientation
            eo = self.edge_orient[i]
            
            # Loop through dofs in edge nbr. i
            for dof in self.edof[i]:
                
                #dof = int(dof)
                # Diagonal 
                A[dof,dof] = -3
                
                # Dirichlet condition (No A modification)
                if self.edge_type[i] == 'd':
                    # Condition is known and constant (constructs basic f)
                    if self.fetch[i] == None:
                        fv = (-1) * self.edge_init[i] / (self.dx**2)
                        if eo == 't':
                            f[dof+self.x_dof] += fv
                            f[dof-1] += fv
                            f[dof+1] += fv
                        elif eo == 'b':
                            f[dof-self.x_dof] += fv
                            f[dof-1] += fv
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
                    if eo == 'l':
                        A[dof, dof-1] = 0
                    elif eo == 'r':
                        if dof != self.nbr_dof-1:
                            A[dof, dof+1] = 0   
                    
                    # Condition is known and constant
                    if self.fetch[i] == None:  
                        f[dof] += (-1) * self.edge_init[i] / self.dx
        
        vertices_dof = np.array([0, self.x_dof-1, 
                                  self.nbr_dof-self.x_dof, self.nbr_dof-1])
        for dof in vertices_dof:
            A[dof,dof] = -2
            
        A = A / (self.dx**2)
        
        # Add known bc:s to v-vector
        v = np.zeros(self.nbr_dof)
        for i in range(len(self.edof)):
            if self.fetch[i] == None and self.edge_type[i] == 'd':
                v[self.edof[i]] = self.edge_init[i]
        
        A_red = A[self.red_dofs][:,self.red_dofs]

        print(f/4)
        
        omega = 0.8
        # Solving loop
        nbr_iter = 10
        for i in range(nbr_iter):
            # Recieve fetched conditions and update f
            f_new = f.copy()
            # Loop through edges
            for j in range(len(self.edof)):
                # Edge nbr. j
                
                # Edge that fetches boundary condition from other process
                if self.fetch[j] != None:
                    cond = comm.recv(source = self.fetch[j])
                    #cond = np.ones(len(self.edof[j])) # Test value
                    
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
                            f_new[self.edof[j]+1] -= fv
                           
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
    dx = 0.5
    # Define region
    #points = np.array([(0,2), (1,2), (1,1), (1,0), (0,0), (0,1)])
    points = np.array([(0,1), (1,1), (1,0), (0,0)])
    edge_type = np.array(['d', 'n', 'd', 'd'])
    fetch = np.array([None, 3, None, None])
    edge_init = np.array([40, None, 5, 15])

    # Create region
    r = region(points, edge_type, fetch, edge_init, dx)
    
    # Solve
    v = r.solve(1)