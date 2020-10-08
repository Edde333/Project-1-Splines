import numpy as np
from mpi4py import MPI
from region import region
import matplotlib.pyplot as plt

if __name__ == "__main__":

    dx = 0.5
    guess = 20

    # Get communicator
    comm = MPI.COMM_WORLD

    # Region 1 (starts iteration)
    if comm.Get_rank() == 1:
        # Define region
        points = np.array([(0,2), (1,2), (1,1), (1,0), (0,0), (0,1)])
        edge_type = np.array(['d', 'd', 'd', 'd', 'd', 'd'])
        fetch = np.array([None, 3, None, None, None, None])
        edge_init = np.array([40, None, 15, 5, None, 15])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        
        # Solve system
        r.solve(comm)

    # Region 2 (bottom left)
    if comm.Get_rank() == 2:
        # Define region
        points = np.array([(0,1), (1,1), (1,0), (0,0)])
        edge_type = np.array(['d', 'n', 'd', 'd'])
        fetch = np.array([None, 1, None, None])
        edge_init = np.array([15, None, 15, 40])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        
        # Solve system
        comm.send(guess, dest=1)
        r.solve(comm)

    # Region 3 (top right)
    if comm.Get_rank() == 3:
        # Define region
        points = np.array([(0,1), (1,1), (1,0), (0,0)])
        edge_type = np.array(['d', 'd', 'd', 'n'])
        fetch = np.array([None, None, None, 1])
        edge_init = np.array([15, 40, 15, None])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        
        # Solve system
        comm.send(guess, dest=1)
        r.solve(comm)

    # Summarize result and plot
    if comm.Get_rank() == 0:
        v1 = comm.recv(source = 1)
        v2 = comm.recv(source = 2)
        v3 = comm.recv(source = 3)
        
        
        
        # Plotting ???