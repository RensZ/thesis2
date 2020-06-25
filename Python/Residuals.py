"""
Created by Rens van der Zwaard on 2020-4-7

Purpose: to plot the residuals of the individual observations for inspection

"""


def f(dir_output, dir_plots, body, no_arcs, useRSW):

    import numpy as np
    import matplotlib.pyplot as plt
    from ToolKit import format_spines
    from os import path

    useAbs = True

    if path.exists(dir_output + "EstimationConsiderCorrelations.dat"):
        c = 2
    else:
        c = 1

    for consider in range(c):

        if consider == 0:
            filestring = ""
            savestring = ""
        else:
            filestring = "Consider"
            savestring = "_consider"

        #residual history
        t = np.genfromtxt(dir_output+"ObservationTimes.dat")
        r = np.genfromtxt(dir_output + "ResidualHistory.dat")

        r_rms = np.sqrt( np.sum(r**2, axis=0) / len(r) )
        if r.ndim == 2:
            i_best = np.where(r_rms == np.min(r_rms))[0]
            r_best = r[:,i_best]
        else:
            r_best = r

        if useRSW:
            filename = "propagatedRSWErrorUsing"+filestring+"CovMatrix.dat"
            legend = ["x", "y", "z"]
            saveRSW = ""
        else:
            filename = "propagatedErrorUsing"+filestring+"CovMatrix.dat"
            legend = ["r", "at", "ct"]
            saveRSW = "_RSW"

        #propagated errors
        if useAbs:
            e_data = np.abs(np.genfromtxt(dir_output + filename))
            yscale = "log"
        else:
            e_data = np.genfromtxt(dir_output + filename)
            yscale = "symlog"
        t_e = e_data[:,0]

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
            gaps_r = np.concatenate([[0],np.where(dt_r>24.0*60.0*60.0)[0],[-1]])
            gaps_e = np.concatenate([[0],np.where(dt_e>24.0*60.0*60.0)[0],[-1]])


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

            if no_arcs > 1:
                t_arc = t_sorted[start_r:end_r]-t_sorted[0]
                r_arc = np.abs(r_sorted[start_r:end_r])
                e_arc = e_data[start_e:end_e,:]
                te_arc = t_e[start_e:end_e]-t_e[0]
            else:
                t_arc = t_sorted-t_sorted[0]
                r_arc = np.abs(r_sorted)
                e_arc = e_data
                te_arc = t_e - t_e[0]

            t_av_r.append(t_sorted[start_r])
            t_av_e.append(t_e[start_e])
            av_r.append(np.mean(r_arc))
            std_r.append(np.std(r_arc))
            av_e[i-1,:] = np.mean(e_arc[:,1:7],axis=0)
            std_e[i-1,:] = np.std(e_arc[:,1:7],axis=0)

            ax = fig.add_subplot(3, no_arcs, i)
            ax.plot(t_arc,r_arc,"ro",markersize=1)
            ax.set_xlabel("t [s]",horizontalalignment='left',x=0.01)
            if i == 1:
                ax.set_ylabel("residual [m]")
            ax.set_yscale("log")
            if no_arcs > 1:
                ax.set_ylim(np.min(np.abs(r_sorted)),np.max(np.abs(r_sorted)))
            format_spines(ax, i, no_arcs)

            ax = fig.add_subplot(3, no_arcs, i+no_arcs)
            for j in range(1, 4):
                ax.plot(te_arc,e_arc[:,j],linewidth=0.75)
            ax.set_xlabel("t [s]", horizontalalignment='left',x=0.01)
            if i == 1:
                ax.set_ylabel("propagated position error [m]")
                ax.legend(legend, loc='upper left')
            ax.set_yscale(yscale)
            plotdata = e_data[:,1:4]
            y_min = np.min(plotdata[np.nonzero(plotdata)])
            y_max = np.max(plotdata)
            ax.set_ylim((y_min, y_max))
            format_spines(ax, i, no_arcs)

            ax = fig.add_subplot(3, no_arcs, i + 2*no_arcs)
            for j in range(4, 7):
                ax.plot(te_arc,e_arc[:,j], linewidth=0.75)
            ax.set_xlabel("t [s]",horizontalalignment='left',x=0.01)
            if i == 1:
                ax.set_ylabel("propagated velocity error [m/s]")
            ax.set_yscale(yscale)
            plotdata = e_data[:,4:7]
            y_min = np.min(plotdata[np.nonzero(plotdata)])
            y_max = np.max(plotdata)
            if not (np.isnan(y_min) or np.isnan(y_max)):
                ax.set_ylim((y_min, y_max))
                format_spines(ax, i, no_arcs)


        plt.tight_layout()
        plt.savefig(dir_plots+body+"_observation_residuals_perarc"+saveRSW+savestring)


        if no_arcs > 1:

            t_av_r = np.asarray(t_av_r)
            t_av_e = np.asarray(t_av_e)

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
            plt.legend(legend, loc='upper left')

            plt.subplot(3,1,3)
            for j in range(3, 6):
                plt.errorbar(t_av_e,av_e[:,j],std_e[:,j],fmt='--o')
            plt.xlabel("time since J2000 [s]")
            plt.ylabel("propagated velocity error [m/s]")
            plt.yscale(yscale)

            plt.tight_layout()
            plt.savefig(dir_plots+body+"_residuals_averages"+saveRSW+savestring)


        #save averages to use as input for thesis_v1.cpp
        if not useRSW:
            save_array = np.vstack((t_av_e, av_e[:, 0], av_e[:, 1], av_e[:, 2], av_e[:, 3], av_e[:, 4], av_e[:, 5])).T
            averages_filename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/error_inputs_"+body+savestring+".txt"
            np.savetxt(averages_filename, save_array, delimiter=",")

