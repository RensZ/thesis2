"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to make a heat map of correlation coefficients between parameters

"""


def f(dir_output, dir_plots, parameters, no_arcs):

    import seaborn as sns
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from os import path

    if path.exists(dir_output + "EstimationConsiderCorrelations.dat"):
        c = 2
    else:
        c = 1

    for consider in range(c):

        if consider == 0:
            correlationmatrix = np.genfromtxt(dir_output+"EstimationCorrelations.dat")
            savestring = ""
        else:
            correlationmatrix = np.genfromtxt(dir_output+"EstimationConsiderCorrelations.dat")
            savestring = "_consider"

        parameters_labels = ['X_Mer', 'Y_Mer', 'Z_Mer', 'Vx_Mer', 'Vy_Mer', 'Vz_Mer']
        for p in parameters[6:]:
            if p == "gamma":
                parameters_labels.append(r"$\gamma$")
            if p == "beta":
                parameters_labels.append(r"$\beta$")
            if p == "TVGP":
                parameters_labels.append(r"$\dot{Gm_{\odot}}/Gm_{\odot}$")
            if p == "J2_Sun":
                parameters_labels.append(r"$J_{2\odot}$")
            if p == "J2_A":
                parameters_labels.append(r"$A_{J2\odot}$")

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
            df = pd.DataFrame(correlationmatrix, index=parameters_labels, columns=parameters_labels)
            df_onlyparameters = pd.DataFrame(correlationmatrix[6:,6:], index=parameters_labels[6:], columns=parameters_labels[6:])
        else:
            df = pd.DataFrame(np.abs(correlationmatrix))
            df_onlyparameters = pd.DataFrame(np.abs(correlationmatrix[6:,6:]))

        sns.set(style="white")
        cmap = sns.diverging_palette(220, 10, as_cmap=True)

        fig = plt.figure(figsize=(10,10))
        if parameters:
            sns.heatmap(df, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                    linewidths=.5, annot=True, cbar=False, fmt='.2f')
        else:
            sns.heatmap(df, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                    linewidths=.1, cbar=False)

        plt.savefig(dir_plots + 'parameter_correlation'+savestring+'.png')


        fig2 = plt.figure(figsize=(4.5,4.5))
        if parameters:
            sns.heatmap(df_onlyparameters, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                    linewidths=.5, annot=True, cbar=False, fmt='.2f')
        else:
            sns.heatmap(df_onlyparameters, cmap=cmap, center=0, vmin=-1, vmax=1, square=True,
                    linewidths=.1, cbar=False)
        plt.tight_layout()
        plt.savefig(dir_plots + 'parameter_correlation_onlyparameters' + savestring + '.png')

    return