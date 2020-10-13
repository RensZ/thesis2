
import numpy as np
import datetime
import matplotlib.pyplot as plt

dir_plots = '/home/rens/Documents/PostProcessing_plots/'

# cycle 24 data
sunspotData = np.genfromtxt("SN_ms_tot_V2.0.txt")

dates24 = []
for i in range(len(sunspotData)):
    dates24.append(datetime.datetime(int(sunspotData[i,0]),int(sunspotData[i,1]), 15))
ssn24 = sunspotData[:,3]

# cycle 25 data
import json
with open("ssn_prediction.json") as f:
    json_input = json.load(f)

dates25 = []
ssn25 = []
for json_string in json_input:
    date_string = json_string["time-tag"]
    year = int(date_string[0:4])
    month = int(date_string[5:])
    if year<2033 and (year>2015):
        dates25.append(datetime.datetime(year,month,15))
        ssn25.append(json_string['predicted_ssn'])

# mission data
MessStart = datetime.datetime(2011,3,18)
MessEnd = datetime.datetime(2015,4,30)
BepiStart = datetime.datetime(2025,12,2)
BepiEnd = datetime.datetime(2027,5,1)


# lets plot
fig = plt.figure(figsize=(10,3))
plt.plot(dates24, ssn24)
plt.plot(dates25, ssn25)
plt.legend(["Cycle 24", "Cycle 25 prediction"])

plt.hlines(np.min(ssn24)-0, MessStart, MessEnd, colors='red', linewidth=5)
plt.hlines(np.min(ssn24)-0, BepiStart, BepiEnd, colors='red', linewidth=5)
plt.text(MessStart-datetime.timedelta(1.0*30.0),np.min(ssn24)+3,"MESSENGER in orbit",color="red")
plt.text(BepiStart-datetime.timedelta(15.0*30.0),np.min(ssn24)+3,"BepiColombo in orbit",color="red")

plt.ylabel("monthly smoothed sunspot number")
plt.ylim(bottom=np.min(ssn24)-1)
plt.grid()
plt.tight_layout()
plt.savefig(dir_plots + "solarcycle.png")
plt.close()
