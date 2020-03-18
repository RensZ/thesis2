"""
Created by Rens van der Zwaard on 2020-3-2

Purpose: to plot the propagation of the bodies to check whether it went alright

"""


def f(dir_output, dir_plots, parameters, bodies):

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    AU = 149597870700.0 #m

    no_bodies = len(bodies)

    for i in range(0, no_bodies):

        body = bodies[i]

        filename = "StatePropagationHistory" + body + ".dat"
        data = np.genfromtxt(dir_output+filename,delimiter=',')
        x = data[:,1]/AU
        y = data[:,2]/AU
        z = data[:,3]/AU

        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')
        ax.plot(x,y,z,linewidth=0.75)
        # ax.plot(x[0:10000],y[0:10000],z[0:10000],linewidth=0.75)
        ax.set_xlabel('x [AU]')
        ax.set_ylabel('y [AU]')
        ax.set_zlabel('z [AU]')

        plt.savefig(dir_plots + body + '_state_history.png')

