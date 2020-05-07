"""
Created by Rens van der Zwaard on 2020-4-23

Purpose: to plot the numerical integration error after integrating backwards

"""


def f(dir_output, dir_plots, body, no_arcs):

    import numpy as np
    import matplotlib.pyplot as plt
    from ToolKit import format_spines

    filename = "StatePropagationHistory"+body+".dat"
    filename_b = "StatePropagationHistory"+body+"Backwards.dat"

    data = np.genfromtxt(dir_output+filename,delimiter=',')
    data_b = np.genfromtxt(dir_output+filename_b, delimiter=',')

    t = data[:,0]
    t_b = data_b[:,0]

    if no_arcs == 1:
        gaps = [0, -1]
        gaps_b = [0, -1]
    else:
        dt = t[1:-1] - t[0:-2]
        gaps = np.concatenate([[0], np.where(dt > 24.0*60.0*60.0)[0], [-1]])
        dt_b = t_b[1:-1] - t_b[0:-2]
        gaps_b = np.concatenate([[0], np.where(dt_b > 24.0*60.0*60.0)[0], [-1]])


    if np.array_equal(t, t_b):
        data_b_interp = data_b[:,1:7]
        print("  timestamps equal for forward and backward integration")
    else:
        for j in range(1, no_arcs + 1):

            start = gaps[j-1]
            start_b = gaps_b[j-1]
            end = gaps[j]
            end_b = gaps_b[j]

            if end == -1:
                t_arc = t[start:]
                t_arc_b = t_b[start_b:]
                data_arc = data[start:, 1:7]
            else:
                t_arc = t[start:end]
                t_arc_b = t_b[start_b:end_b]
                data_arc = data[start:end, 1:7]

            data_arc_b = np.zeros((len(data_arc),6))
            for k in range(0, 6):
                # #linear interpolator
                # data_k = np.interp(t_arc, t_b, data_b[:, k+1])

                # cubic spline interpolator
                from scipy.interpolate import CubicSpline
                if end == -1:
                    cs = CubicSpline(t_arc_b, data_b[start_b:, k+1])
                else:
                    cs = CubicSpline(t_arc_b, data_b[start_b:end_b, k+1])
                data_k = cs(t_arc)

                data_arc_b[:,k] = data_k

            # fig3 = plt.figure(figsize=(16, 10))
            # plt.plot(t_arc_b,data_b[start_b:end_b,1])
            # plt.plot(t_arc, data_arc_b[:,0])
            # plt.savefig(dir_plots + body + str(j) + '_check_interpolation.png')

            if j == 1:
                data_b_interp = data_arc_b
            else:
                data_b_interp = np.vstack([data_b_interp,data_arc_b])


    fig = plt.figure(figsize=(16, 10))
    allerrors = np.abs(data_b_interp - data[:,1:7])
    allerrors_pos = allerrors[:,0:3]
    allerrors_vel = allerrors[:,3:6]

    y_min_pos = np.min(allerrors_pos[np.nonzero(allerrors_pos)])
    y_max_pos = np.max(allerrors_pos)
    y_min_vel = np.min(allerrors_vel[np.nonzero(allerrors_vel)])
    y_max_vel = np.max(allerrors_vel)

    for j in range(1, no_arcs + 1):

        start = gaps[j-1] + 1
        end = gaps[j]
        t_arc = t[start:end]
        error = allerrors[start:end,:]

        ax = fig.add_subplot(2, no_arcs, j)
        for k in range(0, 3):
            ax.plot(t_arc - t[0], error[:,k], linewidth=0.75)
        ax.set_xlabel("t [h]", horizontalalignment='left', x=0.01)
        if j == 1:
            ax.set_ylabel("absolute position integration error [m]")
            ax.legend(["x", "y", "z"])
        ax.set_yscale("log")
        ax.set_ylim(y_min_pos, y_max_pos)
        format_spines(ax, j, no_arcs)

        ax = fig.add_subplot(2, no_arcs, j + no_arcs)
        for k in range(3, 6):
            ax.plot(t_arc - t[0], error[:,k], linewidth=0.75)
        ax.set_xlabel("t [h]", horizontalalignment='left', x=0.01)
        if j == 1:
            ax.set_ylabel("absolute velocity integration error [m/s]")
            ax.legend(["x", "y", "z"])
        ax.set_yscale("log")
        ax.set_ylim(y_min_vel, y_max_vel)
        format_spines(ax, j, no_arcs)

    plt.tight_layout()
    plt.savefig(dir_plots + body + '_state_history_integration_error.png')

    fig2 = plt.figure(figsize=(16, 10))

    allerrors_pos_norm = np.linalg.norm(allerrors_pos,axis=1)
    allerrors_vel_norm = np.linalg.norm(allerrors_vel,axis=1)

    allerrors_pos_norm2 = allerrors_pos_norm[0:-1]
    allerrors_vel_norm2 = allerrors_vel_norm[0:-1]

    y_min_pos = np.min(allerrors_pos_norm2[np.nonzero(allerrors_pos_norm2)])
    y_max_pos = np.max(allerrors_pos_norm2)
    y_min_vel = np.min(allerrors_vel_norm2[np.nonzero(allerrors_vel_norm2)])
    y_max_vel = np.max(allerrors_vel_norm2)

    print("  maximum integration error, position norm [m]:", y_max_pos)
    print("  maximum integration error, velocity norm [m/s]:", y_max_vel)

    for j in range(1, no_arcs + 1):

        start = gaps[j-1]+1
        end = gaps[j]
        t_arc = t[start:end]
        pos_error = allerrors_pos_norm[start:end]
        vel_error = allerrors_vel_norm[start:end]

        ax = fig2.add_subplot(2, no_arcs, j)
        ax.plot(t_arc - t[0], pos_error, linewidth=0.75)
        ax.set_xlabel("t [h]", horizontalalignment='left', x=0.01)
        if j == 1:
            ax.set_ylabel("absolute position integration error [m]")
        ax.set_yscale("log")
        ax.set_ylim(y_min_pos, y_max_pos)
        format_spines(ax, j, no_arcs)

        ax = fig2.add_subplot(2, no_arcs, j + no_arcs)
        ax.plot(t_arc - t[0], vel_error, linewidth=0.75)
        ax.set_xlabel("t [h]", horizontalalignment='left', x=0.01)
        if j == 1:
            ax.set_ylabel("absolute velocity integration error [m/s]")
        ax.set_yscale("log")
        ax.set_ylim(y_min_vel, y_max_vel)
        format_spines(ax, j, no_arcs)

    plt.tight_layout()
    plt.savefig(dir_plots + body + '_state_history_integration_error_norm.png')