"""
Created by Rens van der Zwaard on 2020-3-12

Purpose: to plot each seperate acceleration term

"""


def f(dir_output, dir_plots, dependent_variables):

    import numpy as np
    import matplotlib.pyplot as plt
    import datetime

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

        for i in range(0, no_variables):

            var = dependent_variables[i]

            if var != "exclude":
                j = i*3+1
                acceleration = np.linalg.norm(data[:,j:j+3], axis=1)
                plt.plot(date,acceleration)
                legend.append(var)

        plt.grid(which='major')
        plt.legend(legend, loc='upper right')
        if f == 0:
            plt.savefig(dir_plots + 'dependent_variable_history_reality.png')
        else:
            plt.savefig(dir_plots + 'dependent_variable_history_estimation.png')