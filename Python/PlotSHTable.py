import numpy as np
import matplotlib.pyplot as plt
import datetime

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_SHtable = dir_application + 'Output/SHtable/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'


from os import scandir, mkdir, path
files = [f.path for f in scandir(dir_SHtable)]

baseData = np.genfromtxt(dir_SHtable + "CosCoeffTable_Amp2.515576E-8_Per3.471336E8_Ph-6.686172.dat", delimiter=",")

date = []
for time in baseData[:, 0]:
    date.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=time))

fig = plt.figure(figsize=(12,6))
legend = []
for f in files:
    if f is not dir_SHtable+"CosCoeffTable_Amp2.515576E-8_Per3.471336E8_Ph-6.686172":

        data = np.genfromtxt(f, delimiter=",")[:,1]
        plt.plot(date,data-baseData[:,1],linewidth=0.75)
        legend.append(f[-50:])

plt.legend(legend)
plt.grid()
plt.tight_layout()
plt.savefig(dir_plots+"SHtable.png")