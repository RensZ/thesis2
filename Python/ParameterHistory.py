"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to plot the parameter history along the iterations

"""


def f(dir_output, parameters, dir_plots, no_bodies):

    import numpy as np
    import matplotlib.pyplot as plt

    no_parameters = (len(parameters) - no_bodies * 6)
    subplots = 2*no_bodies + no_parameters

    data = np.genfromtxt(dir_output + "ParameterHistory.dat")

    fig = plt.figure(figsize=(10,10))

    # State history, only look at position and velocity norm
    for i in range(0, no_bodies):
        j = max(i*6-1,0)
        pos = data[j:j+3]
        vel = data[j+3:j+6]
        posnorm = np.linalg.norm(pos,axis=0)
        velnorm = np.linalg.norm(vel,axis=0)

        k = 1+i*2
        plt.subplot(subplots,1,k)
        plt.plot(posnorm)
        plt.yscale('log')
        plt.ylabel('norm(r) body'+ str(i+1))

        plt.subplot(subplots,1,k+1)
        plt.plot(velnorm)
        plt.yscale('log')
        plt.ylabel('norm(V) body'+ str(i+1))

    # Estimatable parameter history
    for i in range(6*no_bodies,len(parameters)):
        par = data[i]
        k = no_bodies*2+1 + i-6*no_bodies
        plt.subplot(subplots,1,k)
        plt.plot(par)
        plt.yscale('log')
        plt.ylabel(str(parameters[i]))
        if i == len(parameters)-1:
            plt.xlabel('number of iterations')

    plt.savefig(dir_plots + 'paremeter_history.png')

    return