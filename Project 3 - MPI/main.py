import numpy as np
from mpi4py import MPI
from region import region
import matplotlib.pyplot as plt
from scipy import interpolate

if __name__ == "__main__":

    dx = 0.05
    guess = 20

    # Get communicator
    comm = MPI.COMM_WORLD

    # Region 1 (starts iteration)
    if comm.Get_rank() == 1:
        # Define region
        points = np.array([(0, 2), (1, 2), (1, 1), (1, 0), (0, 0), (0, 1)])
        edge_type = np.array(['d', 'd', 'd', 'd', 'd', 'd'])
        fetch = np.array([None, 3, None, None, 2, None])
        edge_init = np.array([40, None, 15, 5, None, 15])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        # Solve system
        res = r.solve(comm)
        comm.send(res, dest=0)

    # Region 2 (bottom left)
    if comm.Get_rank() == 2:
        # Define region
        points = np.array([(0, 1), (1, 1), (1, 0), (0, 0)])
        edge_type = np.array(['d', 'n', 'd', 'd'])
        fetch = np.array([None, 1, None, None])
        edge_init = np.array([15, None, 15, 40])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        
        # Solve system
        comm.send(guess, dest=1)
        res = r.solve(comm)
        comm.send(res, dest=0)

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
        res = r.solve(comm)
        comm.send(res, dest=0)

    # Summarize result and plot
    if comm.Get_rank() == 0:
        # Expecting matrix out
        v1 = comm.recv(source=1)
        v2 = comm.recv(source=2)
        v3 = comm.recv(source=3)

        # ----- Plotting ----- #
        prec = 1e-2
        room1mgrid = np.linspace(0, 1, np.shape(v1)[1])
        room2mgrid = np.linspace(0, 2, np.shape(v1)[0])
        # Interpolate between points to get a nice grid:
        v1_int_f = interpolate.interp2d(room1mgrid, room2mgrid, v1)
        v2_int_f = interpolate.interp2d(room1mgrid, room1mgrid, v2)
        v3_int_f = interpolate.interp2d(room1mgrid, room1mgrid, v3)
        v1_int_2d = v1_int_f(np.arange(0, 1, prec), np.arange(0, 2, prec))
        v2_int_2d = v2_int_f(np.arange(0, 1, prec), np.arange(0, 1, prec))
        v3_int_2d = v3_int_f(np.arange(0, 1, prec), np.arange(0, 1, prec))
        black_room = np.zeros((np.shape(v2_int_2d)[0], np.shape(v2_int_2d)[1]))

        # Concatenate 'rooms' to get the full picture. 
        plot_image = np.concatenate((black_room, v2_int_2d))
        plot_image = np.concatenate((plot_image, v1_int_2d), 1)
        right = np.concatenate((v3_int_2d, black_room))
        plot_image = np.concatenate((plot_image, right), 1)
        fig, ax = plt.subplots()
        im = ax.imshow(plot_image, cmap='jet')
        fig.colorbar(im)
        plt.show()
