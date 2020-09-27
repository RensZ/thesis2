"""
Created by Rens van der Zwaard on 2020-9-1

Purpose: to see which asteroids are the largest contributors to the added formal errors

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def f(dir_output, dir_plots, dir_application, parameterlist):

    formalErrorBeforeAsteroids = np.genfromtxt(dir_output + "ObservationFormalEstimationErrorWithConsiderParameters.dat", delimiter=',')
    formalErrorAfterAsteroids = np.genfromtxt(dir_output + "ObservationFormalEstimationErrorWithConsiderIncludingAsteroidsParameters.dat", delimiter=',')
    asteroidNumbers = np.genfromtxt(dir_application + "Input/mu.txt", dtype=int)[:, 0]
    InpopData = np.genfromtxt(dir_application + "Input/AsteroidsINPOP19a.csv", delimiter=',')

    for p in ["gamma", "beta", "J2_Sun"]:

        print("checking consider covariance increase due to asteroids for parameter: ", p)

        if p not in parameterlist:
            continue
        J2index = np.where(np.asarray(parameterlist) == p)[0][0]

        J2beforeAsteroids = formalErrorBeforeAsteroids[J2index]**2
        J2afterAsteroids = formalErrorAfterAsteroids[J2index]**2

        InpopNumbers = InpopData[:,0]
        InpopUncertainties = InpopData[:,2]
        InpopUncertaintiesMatched = np.zeros(len(asteroidNumbers))

        asteroidFolder = dir_output + "considerCovarianceAsteroids/"

        J2formalErrorIncrease = []
        J2formalErrorIncreasePercentage = []
        sumConsiderCovariance = 0

        for a in asteroidNumbers:

            i = np.where(asteroidNumbers == a)[0][0]
            InpopUncertaintiesMatched[i] = InpopUncertainties[np.where(InpopNumbers == a)[0][0]]

            covarianceCurrentAsteroid = np.genfromtxt(asteroidFolder + "asteroid" + str(a) + ".dat")
            considerCovarianceJ2currentAsteroid = covarianceCurrentAsteroid[J2index, J2index]
            sumConsiderCovariance += considerCovarianceJ2currentAsteroid

            J2formalErrorIncrease.append(considerCovarianceJ2currentAsteroid)
            J2formalErrorIncreasePercentage.append(100.0*considerCovarianceJ2currentAsteroid/J2beforeAsteroids)

        indicesSorted = np.argsort(J2formalErrorIncrease)[::-1]
        asteroidNumbersSorted = asteroidNumbers[indicesSorted]
        J2formalErrorIncreaseSorted = np.asarray(J2formalErrorIncrease)[indicesSorted]
        J2formalErrorIncreasePercentageSorted = np.asarray(J2formalErrorIncreasePercentage)[indicesSorted]
        InpopUncertaintiesSorted = InpopUncertaintiesMatched[indicesSorted]

        df = pd.DataFrame(data={"asteroid": asteroidNumbersSorted,
                                "increase": J2formalErrorIncreaseSorted,
                                "p increase": J2formalErrorIncreasePercentageSorted,
                                "apriori unc.": InpopUncertaintiesSorted})
        print(df)
        df.to_latex(dir_plots + 'asteroid_consider_covariance_contributions_'+p+'.txt', header=True, index=True, float_format="%.2e")

        print("before:", J2beforeAsteroids)
        print("after:", J2afterAsteroids)
        totalIncrease = J2afterAsteroids-J2beforeAsteroids
        print("total increase covariance: ", totalIncrease, ",",
              100*(J2afterAsteroids-J2beforeAsteroids)/J2beforeAsteroids, "%")
        print("check: ", sumConsiderCovariance, ",",
              100 * sumConsiderCovariance / J2beforeAsteroids, "%")

        fig = plt.figure(figsize=(12,8))

        plt.plot(InpopUncertaintiesSorted, J2formalErrorIncreaseSorted, "bo", markersize=3)

        ceresIndex = np.where(asteroidNumbersSorted == 1)[0][0]
        plt.plot(InpopUncertaintiesSorted[ceresIndex], J2formalErrorIncreaseSorted[ceresIndex], "ro", markersize=5)

        pallasIndex = np.where(asteroidNumbersSorted == 2)[0][0]
        plt.plot(InpopUncertaintiesSorted[pallasIndex], J2formalErrorIncreaseSorted[pallasIndex], "go", markersize=5)

        vestaIndex = np.where(asteroidNumbersSorted == 4)[0][0]
        plt.plot(InpopUncertaintiesSorted[vestaIndex], J2formalErrorIncreaseSorted[vestaIndex], "co", markersize=5)

        hygieaIndex = np.where(asteroidNumbersSorted == 10)[0][0]
        plt.plot(InpopUncertaintiesSorted[hygieaIndex], J2formalErrorIncreaseSorted[hygieaIndex], "mo", markersize=5)

        plt.hlines(totalIncrease, np.min(InpopUncertaintiesSorted), np.max(InpopUncertaintiesSorted))

        plt.legend(["Other asteroids","Ceres","Pallas","Vesta","Hygiea","Total increase due to asteroids"])

        plt.xlabel("apriori uncertainty asteroid mass [1E-18 AU3/d2]")
        plt.ylabel("added formal variance J2")
        #plt.xscale("log")
        plt.yscale("log")

        plt.savefig(dir_plots + "asteroid_added_covariance_"+p+".png")

    return