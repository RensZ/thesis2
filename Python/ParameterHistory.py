"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to plot the parameter history along the iterations

"""




def f(dir_output, dir_plots, parameters, bodies, json_input):

    import numpy as np
    import matplotlib.pyplot as plt
    from math import ceil
    from ToolKit import Knm
    import pandas as pd


    #### PART 1: PLOT PARAMETER HISTORY

    no_bodies = len(bodies)
    no_parameters = (len(parameters) - no_bodies * 6)

    subplotcolumns = 2

    if ("gamma" in parameters) and ("beta" in parameters):
        subplotrows = ceil((no_bodies*2 + no_parameters+1)/subplotcolumns)
    else:
        subplotrows = ceil((no_bodies*2 + no_parameters)/subplotcolumns)

    data = np.genfromtxt(dir_output + "ParameterHistory.dat")
    truth = np.genfromtxt(dir_output + "TruthParameters.dat")
    parameterFormalSigmas = np.genfromtxt(dir_output + "ObservationFormalEstimationError.dat")

    fig = plt.figure(figsize=(16,10))
    plt.title("True error vs iterations")

    # State history, only look at position and velocity norm
    for i in range(0, no_bodies):
        j = max(i*6-1,0)
        k = 1 + i * 2

        postruth = np.linalg.norm(truth[j:j+3])
        veltruth = np.linalg.norm(truth[j+3:j+6])

        pos = data[j:j+3]
        vel = data[j+3:j+6]
        posnorm = np.linalg.norm(pos,axis=0) - postruth
        velnorm = np.linalg.norm(vel,axis=0) - veltruth

        plt.subplot(subplotrows,subplotcolumns,k)
        plt.yscale('symlog')
        plt.axhline(y=0.0, color='orange', linewidth=0.75, linestyle='--')
        plt.plot(posnorm)
        plt.ylabel('norm(r) body'+ str(i+1))

        plt.subplot(subplotrows,subplotcolumns,k+1)
        plt.yscale('symlog')
        plt.axhline(y=0.0, color='orange', linewidth=0.75, linestyle='--')
        plt.plot(velnorm)
        plt.ylabel('norm(V) body'+ str(i+1))

    # Estimatable parameter history
    for i in range(6*no_bodies,len(parameters)):
        k = no_bodies * 2 + 1 + i - 6 * no_bodies

        partruth = truth[i]
        par = data[i] - partruth

        plt.subplot(subplotrows,subplotcolumns,k)
        plt.axhline(y=0.0,color='orange',linewidth=0.75, linestyle='--')
        plt.plot(par)
        plt.yscale('symlog')
        plt.ylabel(str(parameters[i]))

        if i >= len(parameters)-subplotcolumns:
            plt.xlabel('number of iterations')

        if parameters[i] == "gamma":
            indexgamma = i
            gamma = data[i]
            gammaFormalError = parameterFormalSigmas[i]
        if parameters[i] == "beta":
            indexbeta = i
            beta = data[i]
            betaFormalError = parameterFormalSigmas[i]

        if parameters[i] == "J2_Sun":
            J2 = data[i,-1]*Knm(2,0)
            print("  unnormalized J2 result: ", J2)
        elif parameters[i] == "J4_Sun":
            J4 = data[i,-1]*Knm(4,0)
            print("  unnormalized J4 result: ", J4)

    #if both gamma and beta are included, also plot nordtvedt parameter using the linear relation
    if ("gamma" in parameters) and ("beta" in parameters):
        nordtvedt = 4*beta-gamma-3
        #print("  nordtvedt parameter from gamma and beta:", nordtvedt)

        k = no_bodies * 2 + 1 + len(parameters) - 6 * no_bodies
        plt.subplot(subplotrows, subplotcolumns, k)
        plt.axhline(y=0.0, color='orange', linewidth=0.75, linestyle='--')
        plt.plot(nordtvedt)
        plt.yscale('symlog')
        plt.ylabel("nv_constraint")

    plt.tight_layout()
    plt.savefig(dir_plots + 'paremeter_history.png')


    #### PART 2: PRINT FORMAL ERRORS

    parameters2 = []
    outputFormalSigmas = []
    paperFormalSigmas = []
    ratioFormalSigmas = []

    for i in range(6 * no_bodies, len(parameters)):
        p = parameters[i]
        parameters2.append(p)

        if p == "J2_Sun":
            K = Knm(2,0)
        else:
            K = 1.0

        outputFormalSigmas.append(K*parameterFormalSigmas[i])
        fs = json_input["formalSigma_" + p]
        paperFormalSigmas.append(fs)
        ratioFormalSigmas.append(K*parameterFormalSigmas[i]/fs)

    if ("gamma" in parameters) and ("beta" in parameters) and json_input["useNordtvedtConstraint"]:
        CovMatrix = np.genfromtxt(dir_output + "UnnormalizedCovariance.dat")
        gammaBetaCovariance = CovMatrix[indexgamma,indexbeta]
        nordtvedtFormalError = np.sqrt((4*betaFormalError)**2
                                       +gammaFormalError**2
                                       -2*4*gammaBetaCovariance)
        fs = json_input["formalSigma_Nordtvedt"]
        outputFormalSigmas.append(nordtvedtFormalError)
        paperFormalSigmas.append(fs)
        ratioFormalSigmas.append(nordtvedtFormalError/fs)
        parameters2.append("NordtvedtEq")

    df = pd.DataFrame(data={"parameter":parameters2,
                            "estimation":outputFormalSigmas,
                            "publication":paperFormalSigmas,
                            "ratio e/p":ratioFormalSigmas})
    print(df)

    return