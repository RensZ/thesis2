"""
Created by Rens van der Zwaard on 2020-3-2

Purpose: to plot the propagation of the bodies to check whether it went alright

"""


def f(dir_output, dir_plots, body, no_arcs):

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    AU = 149597870700.0 #m

    filename = "StatePropagationHistory"+body+".dat"
    data = np.genfromtxt(dir_output+filename,delimiter=',')

    t = data[:,0]
    i_sorted = t.argsort()
    t_sorted = t[i_sorted]
    dt = t_sorted[1:-1]-t_sorted[0:-2]
    gaps = np.concatenate([[0],np.where(dt>60.0*60.0*24.0)[0],[-1]])

    print("  minimum step size in entire simulation:", np.min(dt), "seconds")

    x = data[:,1]/AU
    y = data[:,2]/AU
    z = data[:,3]/AU

    x_sorted = x[i_sorted]
    y_sorted = y[i_sorted]
    z_sorted = z[i_sorted]

    axmin = np.min([x,y,z])
    axmax = np.max([x,y,z])

    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca(projection='3d')

    legend = []
    for j in range(1, no_arcs + 1):
        start = gaps[j-1] + 1
        end = gaps[j]

        # dt_arc = dt[start:end]
        # print("minimum step size for arc", str(j), "is", str(np.min(dt_arc)), "seconds" )

        ax.plot(x[start:end],y[start:end],z[start:end],linewidth=0.75)
        legend.append("arc" + str(j))

    ax.set_xlabel('x [AU]')
    ax.set_ylabel('y [AU]')
    ax.set_zlabel('z [AU]')
    ax.set_xlim([axmin, axmax])
    ax.set_ylim([axmin, axmax])
    ax.set_zlim([axmin, axmax])

    plt.legend(legend, loc='upper left')

    # #plot mercury for reference
    # u = np.linspace(0, 2 * np.pi, 100)
    # r_mercury = 1.63 * (10**-5) #AU
    # x = r_mercury * np.outer(np.cos(u), np.sin(u))
    # y = r_mercury * np.outer(np.sin(u), np.sin(u))
    # z = r_mercury * np.outer(np.ones(np.size(u)), np.cos(u))
    # ax.plot_surface(x, y, z, color='black')

    plt.tight_layout()
    plt.savefig(dir_plots + body + '_state_history.png')

