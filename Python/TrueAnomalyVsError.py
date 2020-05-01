"""
Created by Rens van der Zwaard on 2020-4-30

Purpose: to plot the relation between the position/velocity errors and the true anomaly of the Mercury orbiter

"""


def f(dir_output, dir_plots):

    import numpy as np
    import matplotlib.pyplot as plt
    from math import pi

    # import data
    TA_data = np.genfromtxt(dir_output + "TrueAnomaly.dat")
    errors = np.genfromtxt(dir_output + "propagatedRSWErrorUsingCovMatrix.dat")

    t = TA_data[:,0]
    t2 = errors[:,0]

    TA = TA_data[:,1]*180.0/pi
    pos = errors[:,1:4]
    vel = errors[:,4:7]

    pos_norm = np.linalg.norm(pos, axis=1)
    vel_norm = np.linalg.norm(vel, axis=1)

    # fig = plt.figure(figsize=(16,10))
    #
    # plt.subplot(2,2,1)
    # for i in range(0, 3):
    #     plt.plot(TA,pos[:,i],"o",markersize=0.5)
    # plt.xlabel("true anomaly [deg]")
    # plt.ylabel("position error [m]")
    # plt.legend(["x","y","z"])
    # plt.yscale("log")
    #
    # plt.subplot(2,2,2)
    # for i in range(0, 3):
    #     plt.plot(TA,vel[:,i],"o",markersize=0.5)
    # plt.xlabel("true anomaly [deg]")
    # plt.ylabel("velocity error [m/s]")
    # plt.legend(["vx","vy","vz"])
    #
    # plt.subplot(2,2,3)
    # plt.plot(TA, pos_norm, "o", markersize=1)
    # plt.xlabel("true anomaly [deg]")
    # plt.ylabel("velocity error [m/s]")
    # plt.legend(["norm"])
    #
    # plt.subplot(2,2,4)
    # plt.plot(TA, vel_norm, "o", markersize=1)
    # plt.xlabel("true anomaly [deg]")
    # plt.ylabel("velocity error [m/s]")
    # plt.legend(["norm"])
    #
    # plt.tight_layout()
    # plt.savefig(dir_plots+"ErrorsVSTrueAnomaly")

    fig = plt.figure(figsize=(16,10))

    plt.subplot(2,2,1)
    for i in range(0, 3):
        plt.plot(TA,pos[:,i],"o",markersize=0.5)
    plt.xlabel("true anomaly [deg]")
    plt.ylabel("position error [m]")
    plt.legend(["x","y","z"])
    plt.yscale("log")

    plt.subplot(2,2,2)
    for i in range(0, 3):
        plt.plot(TA,vel[:,i],"o",markersize=0.5)
    plt.xlabel("true anomaly [deg]")
    plt.ylabel("velocity error [m/s]")
    plt.legend(["vx","vy","vz"])
    plt.yscale("log")

    plt.subplot(2,2,3)
    plt.plot(TA, pos_norm, "o", markersize=1)
    plt.xlabel("true anomaly [deg]")
    plt.ylabel("position error [m/s]")
    plt.legend(["norm"])
    plt.yscale("log")

    plt.subplot(2,2,4)
    plt.plot(TA, vel_norm, "o", markersize=1)
    plt.xlabel("true anomaly [deg]")
    plt.ylabel("velocity error [m/s]")
    plt.legend(["norm"])
    plt.yscale("log")

    plt.tight_layout()
    plt.savefig(dir_plots+"ErrorsVSTrueAnomaly_log")
