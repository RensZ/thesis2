"""
Created by Rens van der Zwaard on 2020-4-30

Purpose: to plot the relation between the position/velocity errors and the true anomaly of the Mercury orbiter

"""


def f(dir_output, dir_plots, body, no_arcs, useRSW):

    import numpy as np
    import matplotlib.pyplot as plt
    from math import pi

    useAbs = True

    if useRSW:
        filename = "propagatedRSWErrorUsingCovMatrix.dat"
        legend = ["x", "y", "z"]
        saveRSW = ""
    else:
        filename = "propagatedErrorUsingCovMatrix.dat"
        legend = ["R", "A-T", "C-T"]
        saveRSW = "_RSW"

    # import data
    if useAbs:
        TA_data = np.abs(np.genfromtxt(dir_output + filename))
        yscale = "log"
    else:
        TA_data = np.genfromtxt(dir_output + filename)
        yscale = "symlog"
    TA_data_raw = np.genfromtxt(dir_output + "TrueAnomaly.dat")
    errors_raw = np.genfromtxt(dir_output + filename)

    goodrows = []
    for i in range(0, len(errors_raw)):
        if not np.array_equal(errors_raw[i, 1:7], np.zeros(6)):
            goodrows.append(i)
    print(" ",len(errors_raw)-len(goodrows), "rows are taken out due to being zero, this figure should be equal to:", no_arcs*2)

    errors = errors_raw[goodrows,:]
    TA_data = TA_data_raw[goodrows,:]

    t = TA_data[:,0]
    t2 = errors[:,0]
    if len(t) != len(t2):
        print("!! ERROR: lenghts of time lists not equal !!")

    if no_arcs == 1:
        gaps = [0,-1]
    else:
        dt = t[1:-1] - t[0:-2]
        gaps = np.concatenate([[0], np.where(dt > 24.0*60.0*60.0)[0], [len(t)]])

    arcnumbers = []
    for i in range(0,no_arcs):
        currentarc = np.full(gaps[i+1]-gaps[i], i)
        arcnumbers = np.concatenate((arcnumbers, currentarc))



    TA = TA_data[:,1]*180.0/pi
    pos = errors[:,1:4]
    vel = errors[:,4:7]

    pos_norm = np.linalg.norm(pos, axis=1)
    vel_norm = np.linalg.norm(vel, axis=1)

    fig = plt.figure(figsize=(18,10))

    for i in range(0, 3):
        plt.subplot(4, 2, 2*i+1)
        plt.plot(TA,pos[:,i],"o",markersize=0.25)
        plt.xlabel("true anomaly [deg]")
        plt.ylabel("error " + legend[i]+" [m]")
        plt.yscale(yscale)
        plt.grid()

    for i in range(0, 3):
        plt.subplot(4, 2, 2*i+2)
        plt.plot(TA,vel[:,i],"o",markersize=0.25)
        plt.xlabel("true anomaly [deg]")
        plt.ylabel("error " + legend[i]+" [m/s]")
        plt.yscale(yscale)
        plt.grid()

    plt.subplot(4,2,7)
    plt.plot(TA, pos_norm, "o", markersize=0.25)
    plt.xlabel("true anomaly [deg]")
    plt.ylabel("position error [m/s]")
    plt.yscale(yscale)
    plt.grid()

    plt.subplot(4,2,8)
    plt.plot(TA, vel_norm, "o", markersize=0.25)
    plt.xlabel("true anomaly [deg]")
    plt.ylabel("velocity error [m/s]")
    plt.yscale(yscale)
    plt.grid()

    plt.tight_layout()
    plt.savefig(dir_plots+"ErrorsVSTrueAnomaly"+saveRSW)

    if not useRSW:
        save_array = np.vstack((arcnumbers, t, TA, errors[:,1], errors[:,2], errors[:,3], errors[:,4], errors[:,5], errors[:,6])).T
        averages_filename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/trueanomaly_inputs_"+body+".txt"
        np.savetxt(averages_filename, save_array, delimiter=",")

        if no_arcs > 1:
            gaps_filename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/arcindices_"+body+".txt"
            np.savetxt(gaps_filename, np.asarray(gaps), delimiter=",")
