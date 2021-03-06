"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: to plot the parameter history along the iterations

"""




def f(dir_output, dir_plots, parameters, no_bodies, json_input, useformalsigmas):

    import numpy as np
    import matplotlib.pyplot as plt
    from math import ceil
    from ToolKit import Knm
    import pandas as pd
    from os import path

    if path.exists(dir_output + "ObservationFormalEstimationErrorWithConsiderIncludingAsteroidsParameters.dat"):
        c = 3
    elif path.exists(dir_output + "ObservationFormalEstimationErrorWithConsiderParameters.dat"):
        c = 2
    else:
        c = 1

    originalFormalError = np.zeros(len(parameters))

    for consider in range(c):

        if consider == 0:
            CovMatrix = np.genfromtxt(dir_output + "InitialCovarianceMatrix.dat")
            savestring = ""
        elif consider == 1:
            CovMatrix = np.genfromtxt(dir_output + "ConsiderCovarianceMatrix.dat")
            savestring = "_consider"
        else:
            CovMatrix = np.genfromtxt(dir_output + "ConsiderIncludingAsteroidsCovarianceMatrix.dat")
            savestring = "_considerIncludingAsteroids"

        parameterFormalSigmas = np.sqrt(np.diagonal(CovMatrix))

        #### PART 1: PLOT PARAMETER HISTORY

        no_parameters = (len(parameters) - no_bodies * 6)

        subplotcolumns = 2

        if (("gamma" in parameters) or (json_input["gammaIsAConsiderParameter"])) and ("beta" in parameters):
            subplotrows = ceil((no_bodies*2 + no_parameters+1)/subplotcolumns)
        else:
            subplotrows = ceil((no_bodies*2 + no_parameters)/subplotcolumns)

        data = np.genfromtxt(dir_output + "ParameterHistory.dat")
        truth = np.genfromtxt(dir_output + "TruthParameters.dat")

        fig = plt.figure(figsize=(16,10))
        plt.title("True error vs iterations")

        # State history, only look at position and velocity norm
        for i in range(0, no_bodies):
            j = max(i*6-1,0)
            k = 1 + i * 2

            postruth = np.linalg.norm(truth[j:j+3])
            veltruth = np.linalg.norm(truth[j+3:j+6])

            pos = data[j:j+3,:]
            vel = data[j+3:j+6,:]
            posnorm = np.linalg.norm(pos,axis=0) - postruth
            velnorm = np.linalg.norm(vel,axis=0) - veltruth

            plt.subplot(subplotrows,subplotcolumns,k)
            plt.yscale('log')
            plt.axhline(y=0.0, color='orange', linewidth=0.75, linestyle='--')
            plt.plot(np.abs(posnorm))
            plt.grid()
            plt.ylabel('norm(r) body'+ str(i+1))

            plt.subplot(subplotrows,subplotcolumns,k+1)
            plt.yscale('log')
            plt.axhline(y=0.0, color='orange', linewidth=0.75, linestyle='--')
            plt.plot(np.abs(velnorm))
            plt.grid()
            plt.ylabel('norm(V) body'+ str(i+1))

        # Estimatable parameter history
        for i in range(6*no_bodies,len(parameters)):
            k = no_bodies * 2 + 1 + i - 6 * no_bodies

            partruth = truth[i]
            par = data[i] - partruth

            plt.subplot(subplotrows,subplotcolumns,k)
            plt.axhline(y=0.0,color='orange',linewidth=0.75, linestyle='--')
            plt.plot(np.abs(par))
            plt.yscale('log')
            plt.grid()
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
            if parameters[i] == "alpha1":
                indexalpha1 = i
                alpha1 = data[i]
                alpha1FormalError = parameterFormalSigmas[i]
            if parameters[i] == "alpha2":
                indexalpha2 = i
                alpha2 = data[i]
                alpha2FormalError = parameterFormalSigmas[i]

            if (parameters[i] == "J2_Sun") or (parameters[i] == "J2_A"):
                J2 = data[i,-1]*Knm(2,0)
                print("  unnormalized " + parameters[i] + " result: ", J2)
            elif (parameters[i] == "J4_Sun") or (parameters[i] == "J4_A"):
                J4 = data[i,-1]*Knm(4,0)
                print("  unnormalized " + parameters[i] + " result: ", J4)

        #if both gamma and beta are included, also plot nordtvedt parameter using the linear relation
        if (("gamma" in parameters) or (json_input["gammaIsAConsiderParameter"]))  and ("beta" in parameters):

            if json_input["gammaIsAConsiderParameter"]:
                gamma = np.ones(len(beta))

            if ("alpha1" in parameters) and ("alpha2" in parameters):
                nordtvedt = 4.0*beta-gamma-3.0-alpha1-(2.0/3.0)*alpha2
            else:
                nordtvedt = 4.0*beta-gamma-3.0
            print("  nordtvedt parameter from ppn parameters:", nordtvedt)

            k = no_bodies * 2 + 1 + len(parameters) - 6 * no_bodies
            plt.subplot(subplotrows, subplotcolumns, k)
            plt.axhline(y=0.0, color='orange', linewidth=0.75, linestyle='--')
            plt.plot(np.abs(nordtvedt))
            plt.yscale('log')
            plt.grid()
            plt.ylabel("nv_constraint")

        plt.tight_layout()
        plt.savefig(dir_plots + 'paremeter_history.png')


        #### PART 2: PRINT FORMAL ERRORS

        parameters2 = []
        trueErrors = []
        estimatedValues = []
        aPrioriSigmas = []
        outputFormalSigmas = []
        factorOfImprovement = []
        paperFormalSigmas = []
        ratioFormalSigmas = []

        trueErrorOutput = np.genfromtxt(dir_output + "ObservationTrueEstimationError.dat")

        for i in range(6 * no_bodies, len(parameters)):
            p = parameters[i]
            parameters2.append(p)

            if p == "J2_Sun":
                K = Knm(2,0)
            elif p == "J4_Sun":
                K = Knm(4,0)
            else:
                K = 1.0

            estimatedValues.append(K*(truth[i] + trueErrorOutput[i]))
            trueErrors.append(K*trueErrorOutput[i])

            aps = json_input["sigma_" + p]
            outputFormalSigmas.append(K*parameterFormalSigmas[i])
            aPrioriSigmas.append(aps)
            factorOfImprovement.append(aps/(K*parameterFormalSigmas[i]))

            if useformalsigmas:
                fs = json_input["formalSigma_" + p]
                paperFormalSigmas.append(fs)
                ratioFormalSigmas.append(K*parameterFormalSigmas[i]/fs)

            if parameters[i] == "TVGP":
                TVGPindex = i
                bestResidualsIndex = np.argmin(np.abs(data[TVGPindex, :] - trueErrorOutput[TVGPindex]  ))

        if (("gamma" in parameters) or (json_input["gammaIsAConsiderParameter"])) and ("beta" in parameters) and json_input["useNordtvedtConstraint"]:

            if "gamma" not in parameters:
                gammaFormalError = json_input["sigma_gamma"]
                gammaBetaCovariance = 0.0
            else:
                gammaBetaCovariance = CovMatrix[indexgamma,indexbeta]

            #variance according to: https://stats.stackexchange.com/questions/160230/variance-of-linear-combinations-of-correlated-random-variables

            if ("alpha1" in parameters) and ("alpha2" in parameters):
                print("   alphas are estimated parameters")
                if "gamma" not in parameters:
                    gammaAlpha1Covariance = 0.0
                    gammaAlpha2Covariance = 0.0
                else:
                    gammaAlpha1Covariance = CovMatrix[indexgamma,indexalpha1]
                    gammaAlpha2Covariance = CovMatrix[indexgamma,indexalpha2]

                betaAlpha1Covariance = CovMatrix[indexbeta,indexalpha1]
                betaAlpha2Covariance = CovMatrix[indexbeta,indexalpha2]
                alpha1Alpha2Covariance = CovMatrix[indexalpha1,indexalpha2]

                nordtvedtFormalVariance = gammaFormalError**2 \
                                            + (4.0*betaFormalError)**2 \
                                            + alpha1FormalError**2 \
                                            + ((2.0/3.0)*alpha2FormalError)**2 \
                                            + 2.0*(-1.0)*4.0*gammaBetaCovariance \
                                            + 2.0 * (-1.0) * (-1.0) * gammaAlpha1Covariance \
                                            + 2.0 * (-1.0) * (-2.0 / 3.0) * gammaAlpha2Covariance \
                                            + 2.0*4.0*(-1.0)*betaAlpha1Covariance \
                                            + 2.0*4.0*(-2.0/3.0)*betaAlpha2Covariance \
                                            + 2.0*(-1.0)*(-2.0/3.0)*alpha1Alpha2Covariance

            elif json_input["ppnAlphasAreConsiderParameters"]:
                print("   alphas are considered parameters")

                if json_input["outputSubFolderName"] == "Schettino2015":
                    correlationEstimationWithAlphas = np.genfromtxt(
                        "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/Schettino2015_alphas/EstimationCorrelations.dat")
                elif json_input["outputSubFolderName"] == "Imperi2018_nvtrue":
                    correlationEstimationWithAlphas = np.genfromtxt(
                        "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/Imperi2018_nvtrue_flybys_alphas/EstimationCorrelations.dat")
                elif json_input["outputSubFolderName"] == "Imperi2018_nvfalse":
                    correlationEstimationWithAlphas = np.genfromtxt(
                        "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/Imperi2018_nvfalse_flybys_alphas/EstimationCorrelations.dat")
                else:
                    print("ERROR? ", json_input["outputSubFolderName"])

                alpha1FormalError = json_input["sigma_alpha1"]
                alpha2FormalError = json_input["sigma_alpha2"]

                gammaAlpha1Covariance = 0.0 # correlationEstimationWithAlphas[6,8]*gammaFormalError*alpha1FormalError
                gammaAlpha2Covariance = 0.0 # correlationEstimationWithAlphas[6,9]*gammaFormalError*alpha2FormalError
                betaAlpha1Covariance = 0.0 # correlationEstimationWithAlphas[7,8]*betaFormalError*alpha1FormalError
                betaAlpha2Covariance = 0.0 # correlationEstimationWithAlphas[7,9]*betaFormalError*alpha2FormalError
                alpha1Alpha2Covariance = 0.0 # correlationEstimationWithAlphas[8,9]*alpha1FormalError*alpha2FormalError

                print(gammaFormalError**2, (4.0*betaFormalError)**2, alpha1FormalError**2, ((2.0/3.0)*alpha2FormalError)**2 ,
                      -8.0*gammaBetaCovariance, 2.0*gammaAlpha1Covariance, (4.0/3.0)*gammaAlpha2Covariance,
                      -8.0*betaAlpha1Covariance, (-16.0/3.0)*betaAlpha2Covariance, (4.0/3.0)*alpha1Alpha2Covariance)

                nordtvedtFormalVariance = gammaFormalError**2 \
                                            + (4.0*betaFormalError)**2 \
                                            + alpha1FormalError**2 \
                                            + ((2.0/3.0)*alpha2FormalError)**2 \
                                            + 2.0*(-1.0)*4.0*gammaBetaCovariance \
                                            + 2.0 * (-1.0) * (-1.0) * gammaAlpha1Covariance \
                                            + 2.0 * (-1.0) * (-2.0 / 3.0) * gammaAlpha2Covariance \
                                            + 2.0*4.0*(-1.0)*betaAlpha1Covariance \
                                            + 2.0*4.0*(-2.0/3.0)*betaAlpha2Covariance \
                                            + 2.0*(-1.0)*(-2.0/3.0)*alpha1Alpha2Covariance
            # elif json_input["ppnAlphasAreConsiderParameters"]:
            #     print("   alphas are considered parameters")
            #
            #     covarianceOutputConsiderEstimation = np.genfromtxt(dir_output + "CovarianceOutputConsiderCovarianceEstimation.dat")
            #     gammaAlpha1Covariance = covarianceOutputConsiderEstimation[7,8]
            #     gammaAlpha2Covariance = covarianceOutputConsiderEstimation[7,9]
            #     betaAlpha1Covariance = covarianceOutputConsiderEstimation[6,8]
            #     betaAlpha2Covariance = covarianceOutputConsiderEstimation[6,9]
            #     alpha1Alpha2Covariance = covarianceOutputConsiderEstimation[8,9]
            #
            #     alpha1FormalError = json_input["sigma_alpha1"]
            #     alpha2FormalError = json_input["sigma_alpha2"]
            #
            #     nordtvedtFormalVariance = (4.0*betaFormalError)**2 \
            #                                 +gammaFormalError**2 \
            #                                 +alpha1FormalError**2 \
            #                                 +((2.0/3.0)*alpha2FormalError)**2 \
            #                                 -2.0*4.0*gammaBetaCovariance \
            #                                 -2.0*4.0*betaAlpha1Covariance \
            #                                 -2.0*4.0*(2.0/3.0)*betaAlpha2Covariance \
            #                                 +2.0*gammaAlpha1Covariance \
            #                                 +2.0*(2.0/3.0)*gammaAlpha2Covariance \
            #                                 +2.0*(2.0/3.0)*alpha1Alpha2Covariance

            else:
                print("   alphas are neglected")
                nordtvedtFormalVariance = (4.0*betaFormalError)**2 \
                                          +gammaFormalError**2 \
                                          -2.0*4.0*gammaBetaCovariance

            nordtvedtFormalError = np.sqrt(nordtvedtFormalVariance)

            aps = json_input["sigma_Nordtvedt"]

            outputFormalSigmas.append(nordtvedtFormalError)
            factorOfImprovement.append(aps / nordtvedtFormalError)
            aPrioriSigmas.append(aps)
            parameters2.append("NordtvedtEq")
            estimatedValues.append(nordtvedt[bestResidualsIndex])
            trueErrors.append(nordtvedt[bestResidualsIndex])


            if useformalsigmas:
                fs = json_input["formalSigma_Nordtvedt"]
                paperFormalSigmas.append(fs)
                ratioFormalSigmas.append(nordtvedtFormalError / fs)


        if useformalsigmas:
            if consider > 0:
                originalTrueToFormalRatio = np.asarray(trueToFormalRatio)
                percentageIncreaseFormalErrors = 100.0 * (np.asarray(outputFormalSigmas) - originalFormalError) / originalFormalError
                df = pd.DataFrame(data={"parameter":parameters2,
                                        "apriori er.":aPrioriSigmas,
                                        "true er.":trueErrors,
                                        "formal er.": outputFormalSigmas,
                                        "t/f ratio": trueToFormalRatio,
                                        # "est. value:":estimatedValues,
                                        "f/a improvement":factorOfImprovement,
                                        "paper formal er.":paperFormalSigmas,
                                        "f/p ratio":ratioFormalSigmas})
                #df["incr. due to C.C."].map(lambda x: '{0:.2f}'.format(x))
            else:
                originalFormalError = np.asarray(outputFormalSigmas)
                trueToFormalRatio = np.abs(np.asarray(trueErrors) / np.asarray(outputFormalSigmas))
                df = pd.DataFrame(data={"parameter":parameters2,
                                        "apriori er.":aPrioriSigmas,
                                        "true er.":trueErrors,
                                        "formal er.": outputFormalSigmas,
                                        "t/f ratio":trueToFormalRatio,
                                        # "est. value:":estimatedValues,
                                        "f/a improvement":factorOfImprovement,
                                        "paper formal er.":paperFormalSigmas,
                                        "f/p ratio":ratioFormalSigmas})

            df['f/p ratio'] = df['f/p ratio'].map(lambda x: '{0:.3f}'.format(x))

        else:
            if consider > 0:
                originalTrueToFormalRatio = np.asarray(trueToFormalRatio)
                percentageIncreaseFormalErrors = 100.0 * (np.asarray(outputFormalSigmas) - originalFormalError) / originalFormalError
                df = pd.DataFrame(data={"parameter":parameters2,
                                        "apriori er.":aPrioriSigmas,
                                        "true er.":trueErrors,
                                        "formal er.": outputFormalSigmas,
                                        "t/f ratio": trueToFormalRatio,
                                        # "est. value:":estimatedValues,
                                        "f/a improvement":factorOfImprovement})
                #df["incr. due to C.C."].map(lambda x: '{0:.2f}'.format(x))
            else:
                originalFormalError = np.asarray(outputFormalSigmas)
                trueToFormalRatio = np.abs(np.asarray(trueErrors) / np.asarray(outputFormalSigmas))
                df = pd.DataFrame(data={"parameter":parameters2,
                                        "apriori er.":aPrioriSigmas,
                                        "true er.":trueErrors,
                                        "formal er.": outputFormalSigmas,
                                        "t/f ratio":trueToFormalRatio,
                                        # "est. value:":estimatedValues,
                                        "f/a improvement":factorOfImprovement})

        # if consider>0:
        #     percentageIncreaseFormalErrors = 100.0*(np.asarray(outputFormalSigmas)-originalFormalError)/originalFormalError
        #     df["incr. due to C.C."] = percentageIncreaseFormalErrors
        #     df["incr. due to C.C."].map(lambda x: '{0:.2f}'.format(x))
        #     df["t/f ratio"] = originalTrueToFormalRatio
        # else:
        #
        #     df["true er."] = trueErrors
        #     df["t/f ratio"] = trueToFormalRatio


        df['f/a improvement'] = df['f/a improvement'].map(lambda x: '{0:.2f}'.format(x))
        df['t/f ratio'] = df['t/f ratio'].map(lambda x: '{0:.2f}'.format(x))

        df.to_latex(dir_plots + 'formal_errors'+savestring+'.txt', header=True, index=True, float_format="%.2e")
        with pd.option_context('display.max_rows', None,
                               'display.max_columns', None,
                               'precision', 2):
            print(df)




    return