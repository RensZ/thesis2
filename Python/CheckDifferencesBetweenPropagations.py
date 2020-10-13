
import numpy as np
import datetime
import matplotlib.pyplot as plt

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_output = dir_application + 'Output/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/runForVariousAmplitudesDelta/'

dir_output_old = dir_output + 'PaperInputs_reality3_estimation3_runForVariousAmplitudesDelta/'
stateHistory_old = np.genfromtxt(dir_output_old + "StatePropagationHistory24.dat", delimiter=',')

time = stateHistory_old[:, 0]
date = []
for t in time:
    date.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=t))

dir_output_new = dir_output + 'PaperInputs_reality3_estimation3_runForVariousAmplitudesDelta/'
amplitudes2 = np.geomspace(1.0e-11, 2.25e-7, endpoint=True, num=24)

maxPosNormList = []

for a in range(len(amplitudes2)):

    stateHistory_new = np.genfromtxt(dir_output_new + "StatePropagationHistory"+str(a)+".dat", delimiter=',')

    stateDifference = stateHistory_new[:, 1:] - stateHistory_old[:, 1:]
    posDifference = stateDifference[:, 0:3]
    velDifference = stateDifference[:, 3:6]

    posNorm = np.linalg.norm(posDifference, axis=1)
    velNorm = np.linalg.norm(velDifference, axis=1)

    maxPosNorm = np.max(posNorm)
    maxPosNormList.append(maxPosNorm)
    print("for amplitude ", amplitudes2[a], " max difference position norm: ", maxPosNorm)

    fig = plt.figure(figsize=(12, 6))

    plt.subplot(2, 1, 1)
    plt.plot(date, posDifference[:, 0])
    plt.plot(date, posDifference[:, 1])
    plt.plot(date, posDifference[:, 2])
    plt.xlabel("time [s]")
    plt.ylabel("position norm difference [m]")
    plt.legend(["x", "y", "z"])

    plt.subplot(2, 1, 2)
    plt.plot(date, posNorm)
    plt.xlabel("time [s]")
    plt.ylabel("maximum position difference [m]")

    plt.tight_layout()
    plt.savefig(dir_plots + "difference" + str(a) + ".png")

    plt.close('all')


fig2 = plt.figure(figsize=(6, 4))
plt.plot(amplitudes2, maxPosNormList, "o", markersize=2.0, linestyle='dashed', linewidth=0.75)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$\Delta$ $A$")
plt.ylabel("position difference after 20 years [m]")

plt.grid()
plt.tight_layout()
plt.savefig(dir_plots + "maxErrorVsAmplitudeDelta.png")

plt.close('all')

