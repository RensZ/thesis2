
import numpy as np
import matplotlib.pyplot as plt
from ToolKit import Knm
import datetime

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

deltas = np.genfromtxt(dir_application + "Input/AmplitudeInputs3.txt")

dir_mid = dir_application + 'Output/PaperInputs_reality3_estimation3_testAmplitudePartial0/'
print(dir_mid)
observationtimes = np.genfromtxt(dir_mid + "saveErrorSamples.dat", delimiter=",")[:, 0]
date = []
for time in observationtimes:
    date.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=time))

partials_mid = np.genfromtxt(dir_mid + "test_unnormalizedPartialDerivatives.dat")
partialsWrtAmp = partials_mid[:,-2]

fig = plt.figure(figsize=(14,10))
plt.subplot(2,1,1)

ratios = []
legend = []

for i in range(int((len(deltas)-1)/2)):
    #
    # if i==0:
    #     continue

    j = 2*i+1

    delta = (2.7671e-08/Knm(2,0))*(deltas[j+1]-deltas[j])
    print("up ", deltas[j+1], "down ", deltas[j], "delta ", delta)

    dir_down = dir_application + 'Output/PaperInputs_reality3_estimation3_testAmplitudePartial'+str(j)+'/'
    dir_up = dir_application + 'Output/PaperInputs_reality3_estimation3_testAmplitudePartial'+str(j+1)+'/'

    obs_down = np.genfromtxt(dir_down + "ObservationMeasurements.dat")
    obs_up = np.genfromtxt(dir_up + "ObservationMeasurements.dat")
    central_difference = (obs_up - obs_down)/delta

    ratios.append(central_difference / partialsWrtAmp)

    # print(central_difference)
    # print(partialsWrtAmp)
    #
    # print(np.average(np.abs(central_difference/partialsWrtAmp)))
    # print(Knm(2,0))

    plt.plot(date, np.abs(central_difference))
    legend.append("c. dif. for delta " + str(deltas[j+1]))

    np.savetxt(dir_application + "Input/CentralDifferenceDelta"+str(int(deltas[j+1]*100.0))+".txt", -1.0*central_difference)

plt.plot(date, np.abs(partialsWrtAmp), linestyle='--', linewidth=0.75)
legend.append("partial derivative")
plt.yscale("symlog")
plt.ylabel("dz / dA")
plt.legend(legend)
plt.grid()

plt.subplot(2,1,2)
for i in range(len(ratios)):
    plt.plot(date, ratios[i], linewidth=0.75)
    print(i, np.average(ratios[i]))
plt.yscale("symlog")
plt.ylabel("ratio c.d. / p.d.")
plt.legend(legend[:-1])
plt.grid()


# partialsWrtJ2 = partials_mid[:,-1]
# relativeToJ2 = partialsWrtAmp / partialsWrtJ2
#
# period = 3.471336E8
# phase = -6.686172
# simpleSine = np.sin(2.0*np.pi*observationtimes/period + phase)
#
# partialsWrtAmpDividedBySine = partialsWrtAmp / simpleSine
#
# plt.subplot(4,1,3)
# plt.plot(date, partialsWrtJ2)
# plt.plot(date, partialsWrtAmpDividedBySine)
# plt.ylim(bottom=-2.0E11, top=2.0E11)
# plt.grid()
#
#
# plt.subplot(4,1,4)
# plt.plot(date, relativeToJ2)
# plt.plot(date, simpleSine)
# plt.ylim(bottom=-1.5, top=1.5)
# plt.grid()

plt.tight_layout()
plt.savefig(dir_plots+"AmplitudePartialTest.png")


