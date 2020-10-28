
import numpy as np
import datetime
import matplotlib.pyplot as plt

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_output = dir_application + 'Output/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

dir_output_old = dir_output + 'PaperInputs_reality1_estimation1/'
stateHistory_old = np.genfromtxt(dir_output_old + "StatePropagationHistoryMercury.dat", delimiter=',')

time = stateHistory_old[:, 0]
date = []
for t in time:
    date.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=t))

dir_output_new = dir_output + 'PaperInputs_reality1_estimation1changeWithFormalError/'
J4Inputs = np.genfromtxt(dir_application + 'Input/changeWithFormalError.txt')

maxPosNormList = []

fig = plt.figure(figsize=(12, 6))
for a in range(len(J4Inputs)):

    stateHistory_new = np.genfromtxt(dir_output_new + "StatePropagationHistory"+str(a)+".dat", delimiter=',')

    stateDifference = stateHistory_new[:, 1:] - stateHistory_old[:, 1:]
    posDifference = stateDifference[:, 0:3]
    velDifference = stateDifference[:, 3:6]

    posNorm = np.linalg.norm(posDifference, axis=1)
    velNorm = np.linalg.norm(velDifference, axis=1)

    maxPosNorm = np.max(posNorm)
    maxPosNormList.append(maxPosNorm)
    print("for J4 = ", J4Inputs[a], " max difference position norm: ", maxPosNorm)

    # plt.subplot(2, 1, 1)
    # plt.plot(date, posDifference[:, 0])
    # plt.plot(date, posDifference[:, 1])
    # plt.plot(date, posDifference[:, 2])
    # plt.xlabel("time [s]")
    # plt.ylabel("position norm difference [m]")
    # plt.legend(["x", "y", "z"])
    #
    # plt.subplot(2, 1, 2)
    plt.plot(date, posNorm)

plt.xlabel("time [s]")
plt.ylabel("maximum position difference [m]")
#plt.yscale("log")

plt.grid()
plt.legend(["gamma", "beta", "eta", "J2", "TVGP"])
plt.tight_layout()
plt.savefig(dir_plots + "changeWithFormalError.png")

plt.close('all')


# fig2 = plt.figure(figsize=(6, 4))
# plt.plot(J4Inputs[1:], maxPosNormList, "o", markersize=2.0, linestyle='dashed', linewidth=0.75)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r"$J_{4\odot}$")
# plt.ylabel("position difference after 20 years [m]")
#
# plt.grid()
# plt.tight_layout()
# plt.savefig(dir_plots + "maxDifferenceVsJ4.png")

plt.close('all')

