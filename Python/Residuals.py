"""
Created by Rens van der Zwaard on 2020-4-7

Purpose: to plot the residuals of the individual observations for inspection

"""


def f(dir_output, dir_plots, body, no_arcs):

    import numpy as np
    import matplotlib.pyplot as plt
    from ToolKit import format_spines

    useAbs = True

    #residual history
    t = np.genfromtxt(dir_output+"ObservationTimes.dat")/(60.0*60.0)
    r = np.genfromtxt(dir_output + "ResidualHistory.dat")

    r_rms = np.sqrt( np.sum(r**2, axis=0) / len(r) )
    if r.ndim == 2:
        i_best = np.where(r_rms == np.min(r_rms))[0]
        r_best = r[:,i_best]
    else:
        r_best = r

    #propagated errors
    if useAbs:
        e_data = np.abs(np.genfromtxt(dir_output + "propagatedErrorUsingCovMatrix.dat"))
        yscale = "log"
    else:
        e_data = np.genfromtxt(dir_output + "propagatedErrorUsingCovMatrix.dat")
        yscale = "symlog"
    t_e = e_data[:,0]/(60.0*60.0)

    #plot per arc the residuals and propagated errors
    if no_arcs == 1:
        gaps_r = [0,-1]
        gaps_e = [0,-1]
        r_sorted = r_best
        t_sorted = t
    else:
        i_sorted = t.argsort()
        t_sorted = t[i_sorted]
        r_sorted = r_best[i_sorted]
        dt_r = t_sorted[1:-1]-t_sorted[0:-2]
        dt_e = t_e[1:-1]-t_e[0:-2]
        gaps_r = np.concatenate([[0],np.where(dt_r>24.0*30.0)[0],[-1]])
        gaps_e = np.concatenate([[0],np.where(dt_e>24.0*30.0)[0],[-1]])

    fig = plt.figure(figsize=(16, 10))

    t_av_r = []
    t_av_e = []
    av_r = []
    std_r = []
    av_e = np.zeros((no_arcs,6))
    std_e = np.zeros((no_arcs,6))


    for i in range(1,no_arcs+1):

        start_r = gaps_r[i-1] + 1
        start_e = gaps_e[i-1] + 1
        end_r = gaps_r[i]
        end_e = gaps_e[i]

        r_arc = np.abs(r_sorted[start_r:end_r])
        arc = e_data[start_e:end_e,1:7]

        t_av_r.append(t_sorted[start_r])
        t_av_e.append(t_e[start_e])
        av_r.append(np.mean(r_arc))
        std_r.append(np.std(r_arc))
        av_e[i-1,:] = np.mean(arc,axis=0)
        std_e[i-1,:] = np.std(arc,axis=0)

        ax = fig.add_subplot(3, no_arcs, i)
        ax.plot(t_sorted[start_r:end_r]-t_sorted[0],r_arc,"ro",markersize=1)
        ax.set_xlabel("t [h]",horizontalalignment='left',x=0.01)
        if i == 1:
            ax.set_ylabel("residual [m]")
        ax.set_yscale("log")
        ax.set_ylim(np.min(np.abs(r_sorted)),np.max(np.abs(r_sorted)))
        format_spines(ax, i, no_arcs)

        ax = fig.add_subplot(3, no_arcs, i+no_arcs)
        for j in range(1, 4):
            ax.plot(t_e[start_e:end_e]-t_e[0],e_data[start_e:end_e,j],linewidth=0.75)
        ax.set_xlabel("t [h]", horizontalalignment='left',x=0.01)
        if i == 1:
            ax.set_ylabel("propagated position error [m]")
            ax.legend(["x", "y", "z"])
        ax.set_yscale(yscale)
        ax.set_ylim(np.min(e_data[:, 1:4]), np.max(e_data[:, 1:4]))
        format_spines(ax, i, no_arcs)

        ax = fig.add_subplot(3, no_arcs, i + 2*no_arcs)
        for j in range(4, 7):
            ax.plot(t_e[start_e:end_e]-t_e[0], e_data[start_e:end_e, j], linewidth=0.75)
        ax.set_xlabel("t [h]",horizontalalignment='left',x=0.01)
        if i == 1:
            ax.set_ylabel("propagated velocity error [m/s]")
        ax.set_yscale(yscale)
        ax.set_ylim(np.min(e_data[:,4:7]), np.max(e_data[:,4:7]))
        format_spines(ax, i, no_arcs)


    plt.tight_layout()
    plt.savefig(dir_plots+body+"_observation_residuals_perarc")

    t_av_r = 60.0 * 60.0 * np.asarray(t_av_r)
    t_av_e = 60.0*60.0*np.asarray(t_av_e)

    # plot average per arc, errobars are the sigmas
    fig2 = plt.figure(figsize=(16, 10))
    plt.subplot(3,1,1)
    plt.errorbar(t_av_r,av_r,std_r,fmt='--o')
    plt.xlabel("time since J2000 [s]")
    plt.ylabel("residual [m]")
    plt.yscale(yscale)


    plt.subplot(3,1,2)
    for j in range(0, 3):
        plt.errorbar(t_av_e,av_e[:,j],std_e[:,j],fmt='--o')
    plt.xlabel("time since J2000 [s]")
    plt.ylabel("propagated position error [m]")
    plt.yscale(yscale)
    plt.legend(["x", "y", "z"])

    plt.subplot(3,1,3)
    for j in range(3, 6):
        plt.errorbar(t_av_e,av_e[:,j],std_e[:,j],fmt='--o')
    plt.xlabel("time since J2000 [s]")
    plt.ylabel("propagated velocity error [m/s]")
    plt.yscale(yscale)
    plt.legend(["x", "y", "z"])

    plt.tight_layout()
    plt.savefig(dir_plots+body+"_residuals_averages")


    #save averages to use as input for thesis_v1.cpp
    save_array = np.vstack((t_av_e, av_e[:,0], av_e[:,1], av_e[:,2], av_e[:,3], av_e[:,4], av_e[:,5])).T
    averages_filename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/error_inputs_"+body+".txt"
    np.savetxt(averages_filename, save_array, delimiter=",")

