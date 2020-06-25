"""
Created by Rens van der Zwaard on 2020-5-7

Purpose: to plot the interpolated errors that were generated based on the results of the mercury orbiter state estimation

"""


def f(dir_output, dir_plots, body, no_arcs, vehicle):

    import numpy as np
    import matplotlib.pyplot as plt
    from ToolKit import format_spines

    filename = "interpolatedErrorMatrix.dat"
    data = np.genfromtxt(dir_output+filename)

    t = data[:,0]
    norm = np.linalg.norm(data[:,1:4],axis=1)

    fig = plt.figure(figsize=(16, 10))
    ylabel = ["x","y","z"]

    for k in range(0, 3):
        plt.subplot(4,1,k+1)
        plt.plot(t, data[:,k+1], linewidth=0.75)
        plt.xlabel("time since J2000 [s]")
        plt.ylabel(ylabel[k] + " satellite error [m]")
        plt.yscale("log")

    plt.subplot(4,1,4)
    plt.plot(t, norm, linewidth=0.75)
    plt.xlabel("time since J2000 [s]")
    plt.ylabel("norm satellite error [m]")
    plt.yscale("log")

    plt.tight_layout()
    plt.savefig(dir_plots + body + "_" + vehicle + '_interpolated_satellite_error.png')

    return