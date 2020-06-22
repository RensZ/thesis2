
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

output_dir = "/home/rens/Documents/PostProcessing_plots/"

data = np.genfromtxt("JnPublications.csv", dtype=None, delimiter=',', skip_header=1)

print(data)

papers = np.asarray([x[0] for x in data]).astype(str)
years = np.asarray([x[1] for x in data]).astype(int)
methods = np.asarray([x[2] for x in data]).astype(str)
J2values = np.asarray([x[3] for x in data]).astype(float)
J2sigmas = np.asarray([x[4] for x in data]).astype(float)
J4values = np.asarray([x[5] for x in data]).astype(float)
J4sigmas = np.asarray([x[6] for x in data]).astype(float)

plt.figure(figsize=(16,10))
plt.errorbar(years, J2values, yerr=J2sigmas)

plt.savefig(output_dir + "J2values")

plt.close()


