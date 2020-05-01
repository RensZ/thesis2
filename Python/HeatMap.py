"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to make a heat map of correlation coefficients between parameters

"""


def f(dir_output, dir_plots, parameters, no_arcs):

    import seaborn as sns
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    correlationmatrix = np.genfromtxt(dir_output+"EstimationCorrelations.dat")
    cn = np.linalg.cond(correlationmatrix)
    print("  condition number full correlation matrix:", "{:e}".format(cn))

    cm_arcs = correlationmatrix[:no_arcs*6-1, :no_arcs*6-1]
    cm_pars = correlationmatrix[no_arcs*6:, no_arcs*6:]
    if len(cm_pars) > 0:
        cn_arcs = np.linalg.cond(cm_arcs)
        print("  condition number arcs only:", "{:e}".format(cn_arcs))
        cn_pars = np.linalg.cond(cm_pars)
        print("  condition number parameters only:", "{:e}".format(cn_pars))

    if parameters:
        df = pd.DataFrame(correlationmatrix, index=parameters, columns=parameters)
    else:
        df = pd.DataFrame(np.abs(correlationmatrix))

    sns.set(style="white")
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    fig = plt.figure(figsize=(10,10))
    if parameters:
        sns.heatmap(df, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                linewidths=.5, annot=True, cbar=False, fmt='.2f')
    else:
        sns.heatmap(df, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                linewidths=.1, cbar=False)

    plt.savefig(dir_plots + 'parameter_correlation.png')

    return