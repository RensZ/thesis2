
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

output_dir = "/home/rens/Documents/PostProcessing_plots/"

data = np.genfromtxt("JnPublications.csv", dtype=None, delimiter=',', skip_header=1)

print(data)

papers   = np.asarray([x[0] for x in data]).astype(str)
years    = np.asarray([x[1] for x in data]).astype(float)
methods  = np.asarray([x[2] for x in data]).astype(str)
J2values = np.asarray([x[3] for x in data]).astype(float)
J2sigmas = 2.0*np.asarray([x[4] for x in data]).astype(float)
J4values = np.asarray([x[5] for x in data]).astype(float)
J4sigmas = np.asarray([x[6] for x in data]).astype(float)

i_physics = np.where(methods == 'Heliophysics')[0]
i_orbits = np.where(methods == 'Planetary orbits')[0]

for i in range(len(years)):
    if papers[i] != 'Mecheri et al 2009':
        while len(np.where((years[i] == years) == True)[0]) > 1:
            years[i] += 0.2


plt.figure(figsize=(8,3))
plt.errorbar(years[i_physics], J2values[i_physics], yerr=J2sigmas[i_physics],
             marker='o', markersize=2.5, linewidth=0, elinewidth=1, capsize=2, capthick=1)
plt.errorbar(years[i_orbits], J2values[i_orbits], yerr=J2sigmas[i_orbits],
             marker='o', markersize=2.5, linewidth=0, elinewidth=1, capsize=2, capthick=1)

# for i, txt in enumerate(papers):
#     plt.annotate(txt, (years[i], J2values[i]))

plt.xlabel('Year of publication')
plt.ylabel('J2 value [-]')
plt.ylim((1.5E-7,3.0E-7))
plt.legend(['Heliophysics','Planetary orbits'], loc='upper left')
plt.grid(axis='y')


plt.tight_layout()
plt.savefig(output_dir + "J2values")

plt.close()


