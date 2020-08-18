"""
Created by Rens van der Zwaard on 2020-3-12

Purpose: to plot each seperate acceleration term, with planets grouped

"""


def f(dir_output, dir_plots, dependent_variables):

    import numpy as np
    import matplotlib.pyplot as plt
    import datetime

    planets = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon"]
    no_variables = len(dependent_variables)

    for f in range(2):

        if f == 0:
            filename = "DependentVariablesHistoryReality.dat"
        else:
            filename = "DependentVariablesHistoryFinalIteration.dat"

        data = np.genfromtxt(dir_output + filename, delimiter=',')
        time = data[:,0]

        date = []
        for t in time:
            date.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=t))

        fig = plt.figure(figsize=(16, 10))

        plt.xlabel('time')
        plt.ylabel('acceleration norm [m/s2]')
        plt.yscale('log')
        legend = []
        planetAcceleration = np.zeros( (len(data),3) )

        for i in range(0, no_variables):

            var = dependent_variables[i]
            j = i * 3 + 1

            if var in planets:
                currentPlanetAcc = data[:, j:j + 3]
                planetAcceleration += currentPlanetAcc
            else:
                if var != "exclude":
                    acceleration = np.linalg.norm(data[:, j:j + 3], axis=1)
                    plt.plot(date,acceleration)
                    legend.append(var)

        plt.plot(date,np.linalg.norm(planetAcceleration, axis=1))
        legend.append("Planets")

        plt.grid(which='major')
        plt.legend(legend, loc='upper right')
        if f == 0:
            plt.savefig(dir_plots + 'dependent_variable_history_planetsgrouped_reality.png')
        else:
            plt.savefig(dir_plots + 'dependent_variable_history_planetsgrouped_estimation.png')