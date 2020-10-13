

import numpy as np

amplitudes = np.geomspace(1.0e-10, 2.25e-7, endpoint=True, num=20)
amplitudes2 = np.geomspace(1.0e-11, 2.25e-7, endpoint=True, num=24)
valuesJ4 = np.geomspace(1.0e-10, 1.0E-5, endpoint=True, num=24)

amplitudesDelta = np.hstack((amplitudes2,[0.0]))
J4Inputs = np.hstack(([0.0],valuesJ4))

print(amplitudes)
print(amplitudesDelta)
print(J4Inputs)

np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/AmplitudeInputs.txt", amplitudes)
np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/AmplitudeInputs2.txt", amplitudesDelta)
np.savetxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/J4Inputs.txt", J4Inputs)