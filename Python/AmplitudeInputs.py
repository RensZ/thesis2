

import numpy as np

amplitudes = np.geomspace(1.0e-10, 2.25e-7, endpoint=True, num=20)
amplitudes2 = np.geomspace(1.0e-11, 2.25e-7, endpoint=True, num=24)
valuesJ4 = np.geomspace(1.0e-10, 1.0E-5, endpoint=True, num=24)

amplitudesDelta = np.hstack((amplitudes2,[0.0]))
J4Inputs = np.hstack(([0.0],valuesJ4))

percentages = [1.0, 5.0, 10.0, 25.0, 50.0]
amplitudesPercentageIncrease = [0.0]
for p in percentages:
    amplitudesPercentageIncrease.append(p/-100.0)
    amplitudesPercentageIncrease.append(p/100.0)

print(amplitudes)
print(amplitudesDelta)
print(amplitudesPercentageIncrease)
print(J4Inputs)

np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/AmplitudeInputs.txt", amplitudes)
np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/AmplitudeInputs2.txt", amplitudesDelta)
np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/AmplitudeInputs3.txt", amplitudesPercentageIncrease)
np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/J4Inputs.txt", J4Inputs)