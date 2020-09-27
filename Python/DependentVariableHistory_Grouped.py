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

    arraysWithJ2Accelerations = []

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
                    if var == "Sun_J2":
                        arraysWithJ2Accelerations.append(acceleration)


        plt.plot(date,np.linalg.norm(planetAcceleration, axis=1))
        legend.append("Planets")

        plt.grid(which='major')
        plt.legend(legend, loc='upper right')
        if f == 0:
            plt.savefig(dir_plots + 'dependent_variable_history_planetsgrouped_reality.png')
        else:
            plt.savefig(dir_plots + 'dependent_variable_history_planetsgrouped_estimation.png')

    # plot difference in J2 acceleration in reality and estimation
    differenceInJ2Acceleration = arraysWithJ2Accelerations[1] - arraysWithJ2Accelerations[0]
    fig2 = plt.figure(figsize=(16, 10))
    plt.plot(date, differenceInJ2Acceleration)
    plt.xlabel('time')
    plt.ylabel('difference J2 acceleration norm [m/s2]')
    # plt.yscale('log')

    observationTimes = np.genfromtxt(dir_output+"interpolatedErrorMatrix.dat")[:,0]
    observationDates = []
    for t in observationTimes:
        observationDates.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=t))
    plt.plot(observationDates, np.mean(differenceInJ2Acceleration)*np.ones(len(observationDates)), "ro", markersize=0.75)

    plt.savefig(dir_plots + 'dependent_variable_J2_difference.png')