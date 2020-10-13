import numpy as np
import matplotlib.pyplot as plt
import datetime

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_SHtable = dir_application + 'Output/SHtable/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'


from os import scandir, mkdir, path
files = [f.path for f in scandir(dir_SHtable)]

fig = plt.figure(figsize=(12,6))
for f in files:

    data = np.genfromtxt(f, delimiter=",")

    date = []
    for time in data[:,0]:
        date.append(datetime.datetime(2000, 1, 1, 12, 0) + datetime.timedelta(seconds=time))

    plt.plot(date,data[:,1],linewidth=0.75)

plt.grid()
plt.tight_layout()
plt.savefig(dir_plots+"SHtable.png")