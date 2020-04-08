"""
Created by Rens van der Zwaard on 2020-4-7

Purpose: to plot the residuals of the individual observations for inspection

"""


def f(dir_output, dir_plots, bodies):

    import numpy as np
    import matplotlib.pyplot as plt

    t_filename = "ObservationTimes.dat"
    t = np.genfromtxt(dir_output+t_filename)

    r_filename = "ResidualHistory.dat"
    r = np.genfromtxt(dir_output+r_filename)

    r_rms = np.sqrt( np.sum(r**2, axis=0) / len(r) )
    i_best = np.where(r_rms == np.min(r_rms))[0]

    r_best = r[:,i_best]

    fig = plt.figure(figsize=(16,10))
    plt.plot(t,r_best,"ro", markersize=1)
    plt.xlabel("time since J2000 [s]")
    plt.ylabel("residual [m]")
    plt.yscale("log")

    plt.savefig(dir_plots+bodies[0]+"_observation_residuals_bestestimate")
