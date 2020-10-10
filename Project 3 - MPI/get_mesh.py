import numpy as np
import sys

def get_mesh(points, dx, fetch):
    
    # Find x- and y-range of region 
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
    y_dof = int(round(height / dx + 1))
    nbr_dof = x_dof * y_dof

    # Find which dof:s belong to the edges
    edof = [np.array([], dtype = int)] * len(points)
    edge_orient = [None] * len(points)
    for i in range(len(edof)):
        # Edge nbr. i-1
        e = np.array([points[i-1], points[i]])

        # Horisontal edge 
        if e[0,1] == e[1,1]:

            # Top edge (i = 0 in v_ij)
            if e[0,1] == y_max:
                edge_orient[i-1] = 't'

                edge_range = np.array([e[0,0]-x_min, e[1,0]-x_min]) / base
                dof_range = (x_dof-1) * edge_range
                edof[i-1] = np.arange(round(dof_range[0]), round(dof_range[1])+1)
                
            # Bottom edge (i = y_dof in v_ij)
            else:
                edge_orient[i-1] = 'b'
                
                edge_range = np.array([e[1,0]-x_min, e[0,0]-x_min]) / base
                dof_range = ((x_dof-1) * edge_range) + (nbr_dof - x_dof)
                edof[i-1] = np.arange(round(dof_range[0]), round(dof_range[1])+1)

        # Vertical edge
        elif e[0,0] == e[1,0]:

            # Left edge (j = 0 in v_ij)
            if e[0,0] == x_min:
                edge_orient[i-1] = 'l'
                
                edge_range = np.array([y_max-e[1,1], y_max-e[0,1]]) / height
                dof_range = ((nbr_dof-x_dof) * edge_range)
                edof[i-1] = np.arange(round(dof_range[0]), round(dof_range[1])+1, x_dof)

            # Right edge    
            else:
                edge_orient[i-1] = 'r'
                edge_range = np.array([y_max-e[0,1], y_max-e[1,1]]) / height
                dof_range = ((nbr_dof-x_dof) * edge_range) + (x_dof-1)
                edof[i-1] = np.arange(round(dof_range[0]), round(dof_range[1])+1, x_dof)
        
        else:
            raise Exception("points vector do not specify a rectangle")
            
    # Handle dofs between edges
    for i in range(len(edof)):
        # If we are fetching we have dof priority
        
        if fetch[i-1] != None and fetch[i] != None:
            raise Exception("Two consecutive edges can't both fetch")
        
        # Edge (i-1) fetches (Edge i deletes common dof)
        if fetch[i-1] != None:
            if edge_orient[i] == 't' or edge_orient[i] == 'r':
                edof[i] = np.delete(edof[i], 0) 
            else: # bottom or right edge
                edof[i] = np.delete(edof[i], -1)
                
        # Edge i fetches or no edge fetches (Edge (i-1) deletes common dof)
        else: 
            if edge_orient[i-1] == 't' or edge_orient[i-1] == 'r':
                edof[i-1] = np.delete(edof[i-1], -1) 
            else: # bottom or right edge
                edof[i-1] = np.delete(edof[i-1], 0)
                
    return x_dof, y_dof, edof, edge_orient

if __name__ == "__main__":
    dx = 0.5
    points = np.array([(0, 2), (1, 2), (1, 1), (1, 0), (0, 0), (0, 1)])
    fetch = np.array([None, 3, None, None, 2, None])
    x_dof, y_dof, edof, edge_orient = get_mesh(points, dx, fetch)
    print(edof)
    
    
    