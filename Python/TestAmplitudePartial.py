
import numpy as np
import matplotlib.pyplot as plt
from ToolKit import Knm

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

delta = (2.7671e-08/Knm(2,0))/10.0

dir_down = dir_application + 'Output/PaperInputs_reality3_estimation3_testAmplitudePartial0/'
dir_mid = dir_application + 'Output/PaperInputs_reality3_estimation3_testAmplitudePartial1/'
dir_up = dir_application + 'Output/PaperInputs_reality3_estimation3_testAmplitudePartial2/'

observationtimes = np.genfromtxt(dir_mid + "saveErrorSamples.dat",delimiter=",")[:,0]

obs_down = np.genfromtxt(dir_down + "ObservationMeasurements.dat")
obs_up = np.genfromtxt(dir_up + "ObservationMeasurements.dat")
central_difference = (obs_up - obs_down)/(2.0*delta)

partials_mid = np.genfromtxt(dir_mid + "test_unnormalizedPartialDerivatives.dat")
partialsWrtAmp = partials_mid[:,-2]

print(central_difference)
print(partialsWrtAmp)

print(np.average(np.abs(central_difference/partialsWrtAmp)))
print(Knm(2,0))

fig = plt.figure(figsize=(14,8))
plt.subplot(2,1,1)
plt.plot(observationtimes, central_difference)
plt.plot(observationtimes, partialsWrtAmp)
plt.legend(["central difference","partial derivative"])
plt.yscale("symlog")
plt.grid()

plt.subplot(2,1,2)
plt.plot(observationtimes, central_difference/partialsWrtAmp)
plt.yscale("symlog")
plt.grid()

plt.tight_layout()
plt.savefig(dir_plots+"AmplitudePartialTest.png")
