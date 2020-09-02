

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("scores.txt", delimiter=" ")
print(data)

x = []
averages = []

for i in range(int(len(data)/10.0)):
    x.append(i*30+15)
    averages.append(np.average(data[i*30:i*30+30]))


print(range(len(data)))
z = np.polyfit(range(len(data)), data, 1)
p = np.poly1d(z)
print("y=%.6fx+(%.6f)"%(z[0],z[1]))


fig = plt.figure()
plt.plot(x,averages)
plt.plot(p(range(len(data))),"r--")
plt.plot([0,len(data)],[1.2,1.2],"g--")

plt.savefig("scores.png")
