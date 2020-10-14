import numpy as np
from mpi4py import MPI
from region import region
import matplotlib.pyplot as plt
from scipy import interpolate
import time as t

if __name__ == "__main__":

    dx = 0.01
    guess = 20
    radiator_heat = 40
    wall_heat = 15
    window_heat = 5

    # Get communicator
    comm = MPI.COMM_WORLD
    
    start_time = t.time()

    # Region 1 (starts iteration)
    if comm.Get_rank() == 1:
        # Define region
        points = np.array([(0, 2), (1, 2), (1, 1), (1, 0), (0, 0), (0, 1)])
        edge_type = np.array(['d', 'd', 'd', 'd', 'd', 'd'])
        fetch = np.array([None, 3, None, None, 2, None])
        edge_init = np.array([radiator_heat, None, wall_heat, window_heat, None, wall_heat])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        print("Process 1 define region:", t.time()-start_time, " s")
        # Solve system
        res = r.solve(comm)
        comm.send(res, dest=0)

    # Region 2 (bottom left)
    if comm.Get_rank() == 2:
        # Define region
        points = np.array([(0, 1), (1, 1), (1, 0), (0, 0)])
        edge_type = np.array(['d', 'n', 'd', 'd'])
        fetch = np.array([None, 1, None, None])
        edge_init = np.array([wall_heat, None, wall_heat, radiator_heat])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        print("Process 2, Define region: ", t.time()-start_time, " s")
        
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
        edge_init = np.array([wall_heat, radiator_heat, wall_heat, None])

        # Create region
        r = region(points, edge_type, fetch, edge_init, dx)
        print("Process 3, Define region: ", t.time()-start_time, " s")
        
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

        t0 = t.time()
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
        
        # Where is the temperature between 20 and 25?
        fig, ax = plt.subplots()
        im = ax.imshow((plot_image <= 25) * (plot_image >= 20), cmap='jet')
        
        plt.show()
        
        t1 = t.time()
        
        print("Plotting:", t1-t0, "s")
        print("Total time consumption:", t1-start_time, "s")
