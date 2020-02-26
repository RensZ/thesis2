"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to make a heat map of correlation coefficients between parameters

"""


def f(dir_output, parameters, dir_plots):

    import seaborn as sns
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    correlationmatrix = np.genfromtxt(dir_output+"EstimationCorrelations.dat")
    df = pd.DataFrame(correlationmatrix, index=parameters, columns=parameters)

    sns.set(style="white")
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    fig = plt.figure(figsize=(10,10))
    sns.heatmap(df, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                linewidths=.5, annot=True, cbar=False, fmt='.2f')
    plt.savefig(dir_plots + 'parameter_correlation.png')

    return