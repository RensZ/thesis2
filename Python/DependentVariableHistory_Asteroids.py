"""
Created by Rens van der Zwaard on 2020-3-12

Purpose: to plot each seperate acceleration term

"""


def f(dir_output, dir_plots):

    import numpy as np
    import matplotlib.pyplot as plt
    import datetime

    planets = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon"]
    first25asteroids = [1, 2, 4, 29, 16, 15, 10, 7, 6, 3, 345, 12, 442, 27, 287, 43, 84, 30, 536, 163, 105, 554, 313, 115, 230]

    dependent_variables = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon"]
    for a in first25asteroids:
        dependent_variables.append("Asteroid"+str(a))
    dependent_variables.append("Sun_CG")
    dependent_variables.append("exclude") #J1
    dependent_variables.append("Sun_J2")
    dependent_variables.append("Sun_SS")
    dependent_variables.append("Sun_SS_α")

    no_variables = len(dependent_variables)

    for j in range(2):

        if j == 0:
            filename = "DependentVariablesHistoryReality.dat"
        else:
            filename = "DependentVariablesHistoryFinalIteration.dat"

        from os.path import isfile
        if not isfile(dir_output + filename):
            print("following file not found, skipping this step: ", filename)
            return

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
        asteroidAcceleration = np.zeros( (len(data),3) )

        for i in range(0, no_variables):

            var = dependent_variables[i]
            j = i * 3 + 1

            if var in planets:
                currentPlanetAcc = data[:, j:j + 3]
                planetAcceleration += currentPlanetAcc
            elif var[0:3] == "Ast":
                currentAsteroidAcc = data[:, j:j + 3]
                asteroidAcceleration += currentAsteroidAcc
            else:
                if var != "exclude":
                    acceleration = np.linalg.norm(data[:, j:j + 3], axis=1)
                    plt.plot(date, acceleration)
                    legend.append(var)

        plt.plot(date, np.linalg.norm(planetAcceleration, axis=1))
        legend.append("Planets")

        plt.plot(date, np.linalg.norm(asteroidAcceleration, axis=1))
        legend.append("Asteroids")

        plt.grid(which='major')
        plt.legend(legend, loc='upper right')
        if j == 0:
            plt.savefig(dir_plots + 'dependent_variable_history_planetsgrouped_reality_asteroids.png')
        else:
            plt.savefig(dir_plots + 'dependent_variable_history_planetsgrouped_estimation_asteroids.png')