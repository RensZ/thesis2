"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to plot the parameter history along the iterations

"""


def f(dir_output, parameters, dir_plots, no_bodies):

    import numpy as np
    import matplotlib.pyplot as plt
    from math import ceil
    from matplotlib.ticker import MaxNLocator

    no_parameters = (len(parameters) - no_bodies * 6)

    subplotcolumns = 2
    subplotrows = ceil( (no_bodies*2 + no_parameters)/subplotcolumns )

    data = np.genfromtxt(dir_output + "ParameterHistory.dat")
    truth = np.genfromtxt(dir_output + "TruthParameters.dat")

    fig = plt.figure(figsize=(16,10))

    # State history, only look at position and velocity norm
    for i in range(0, no_bodies):
        j = max(i*6-1,0)
        k = 1 + i * 2

        pos = data[j:j+3]
        vel = data[j+3:j+6]
        posnorm = np.linalg.norm(pos,axis=0)
        velnorm = np.linalg.norm(vel,axis=0)
        postruth = np.linalg.norm(truth[j:j+3])
        veltruth = np.linalg.norm(truth[j+3:j+6])

        plt.subplot(subplotrows,subplotcolumns,k)
        plt.axhline(y=postruth, color='orange', linewidth=0.75, linestyle='--')
        plt.plot(posnorm)
        plt.yscale('log')
        plt.ylabel('norm(r) body'+ str(i+1))

        plt.subplot(subplotrows,subplotcolumns,k+1)
        plt.axhline(y=veltruth, color='orange', linewidth=0.75, linestyle='--')
        plt.plot(velnorm)
        plt.yscale('log')
        plt.ylabel('norm(V) body'+ str(i+1))

    # Estimatable parameter history
    for i in range(6*no_bodies,len(parameters)):
        par = data[i]
        k = no_bodies*2+1 + i-6*no_bodies
        plt.subplot(subplotrows,subplotcolumns,k)
        plt.axhline(y=truth[i],color='orange',linewidth=0.75, linestyle='--')
        plt.plot(par)
        plt.yscale('log')
        plt.ylabel(str(parameters[i]))
        if i >= len(parameters)-subplotcolumns:
            plt.xlabel('number of iterations')

    plt.tight_layout()
    plt.savefig(dir_plots + 'paremeter_history.png')

    return