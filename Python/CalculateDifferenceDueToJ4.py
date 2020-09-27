
import numpy as np
import matplotlib.pyplot as plt

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_output = dir_application + 'Output/' + "/"
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

ps = "Xu2017"  #"Antia2008" #"Ajab2008"

stateHistory1 = np.genfromtxt(dir_output + ps + "_reality1_estimation1/StatePropagationHistoryMercury.dat",delimiter=',')
stateHistory2 = np.genfromtxt(dir_output + ps + "_reality2_estimation2/StatePropagationHistoryMercury.dat",delimiter=',')

time = stateHistory1[:,0]

stateDifference = stateHistory2[:,1:] - stateHistory1[:,1:]

posDifference = stateDifference[:,0:3]
velDifference = stateDifference[:,3:6]

posNorm = np.linalg.norm(posDifference, axis=1)
velNorm = np.linalg.norm(velDifference, axis=1)

print("max difference position norm: ", np.max(posNorm))
print("max difference velocity norm: ", np.max(velNorm))


fig = plt.figure(figsize=(12,6))

plt.subplot(2,1,1)
plt.plot(time, posDifference[:,0])
plt.plot(time, posDifference[:,1])
plt.plot(time, posDifference[:,2])
plt.xlabel("time [s]")
plt.ylabel("position norm difference [m]")
plt.legend(["x","y","z"])

plt.subplot(2,1,2)
plt.plot(time, posNorm)
plt.xlabel("time [s]")
plt.ylabel("position norm difference [m]")

plt.tight_layout()
plt.savefig(dir_plots + "impact_of_J4_on_position_"+ps+".png")

plt.close('all')