"""
Created by Rens van der Zwaard on 2020-4-23

Purpose: to plot the numerical integration error after integrating backwards

"""


def f(dir_output, dir_plots, body, no_arcs):

    import numpy as np
    import matplotlib.pyplot as plt
    from os import path
    from ToolKit import format_spines

    for m in range(2):

        if m == 0:
            filename_spice = "spiceStatesAtPropagationTimes" + body + ".dat"
            output_string = "_wrt_SPICE"
            legend = ["integration", "SPICE"]
        else:
            filename_spice = "StatePropagationHistory" + body + "Backwards.dat"
            output_string = "_wrt_backwards_integration"
            legend = ["forward", "backward"]


        if not path.exists(dir_output + filename_spice):
            print(dir_output + filename_spice + " not found, returning")
        else:

            filename = "StatePropagationHistory"+body+".dat"

            data = np.genfromtxt(dir_output+filename,delimiter=',')
            spice = np.genfromtxt(dir_output+filename_spice, delimiter=',')

            t = data[:,0]
            t_spice = spice[:,0]

            if not np.array_equal(t, t_spice):
                print("  ERROR: timestamps not equal of integration and spice data! ")
                return


            fig0 = plt.figure(figsize=(16, 10))

            for k in range(0, 6):
                plt.subplot(6,1,k+1)
                plt.plot(t - t[0], data[:,k+1], linewidth=0.75)
                plt.plot(t - t[0], spice[:,k+1], linewidth=0.75)
                plt.xlabel("t [s]", horizontalalignment='left', x=0.01)
                plt.ylabel("state entry "+ str(k+1))
                plt.legend(legend, loc='upper left')

            plt.tight_layout()
            plt.savefig(dir_plots + body + '_comparison_integrated_state'+output_string+'.png')

            allerrors = np.abs(spice[:,1:7] - data[:,1:7])
            allerrors_pos = allerrors[:,0:3]
            allerrors_vel = allerrors[:,3:6]

            y_min_pos = np.min(allerrors_pos[np.nonzero(allerrors_pos)])
            y_max_pos = np.max(allerrors_pos)
            y_min_vel = np.min(allerrors_vel[np.nonzero(allerrors_vel)])
            y_max_vel = np.max(allerrors_vel)


            fig = plt.figure(figsize=(16, 10))

            plt.subplot(2,1,1)
            for k in range(0, 3):
                plt.plot(t - t[0], allerrors_pos[:,k], linewidth=0.75)
            plt.xlabel("t [h]", horizontalalignment='left', x=0.01)
            plt.ylabel("absolute position integration error [m]")
            plt.legend(["x", "y", "z"], loc='upper left')
            plt.yscale("log")
            plt.ylim(y_min_pos, y_max_pos)


            plt.subplot(2,1,2)
            for k in range(0, 3):
                plt.plot(t - t[0], allerrors_vel[:,k], linewidth=0.75)
            plt.xlabel("t [h]", horizontalalignment='left', x=0.01)
            plt.ylabel("absolute velocity integration error [m]")
            plt.legend(["x", "y", "z"], loc='upper left')
            plt.yscale("log")
            plt.ylim(y_min_vel, y_max_vel)

            plt.tight_layout()
            plt.savefig(dir_plots + body + '_state_history_integration_error'+output_string+'.png')



            allerrors_pos_norm = np.linalg.norm(allerrors_pos,axis=1)
            allerrors_vel_norm = np.linalg.norm(allerrors_vel,axis=1)

            y_min_pos = np.min(allerrors_pos_norm[np.nonzero(allerrors_pos_norm)])
            y_max_pos = np.max(allerrors_pos_norm)
            y_min_vel = np.min(allerrors_vel_norm[np.nonzero(allerrors_vel_norm)])
            y_max_vel = np.max(allerrors_vel_norm)

            fig2 = plt.figure(figsize=(16, 10))

            plt.subplot(2,1,1)
            plt.plot(t - t[0], allerrors_pos_norm, linewidth=0.75)
            plt.xlabel("t [h]", horizontalalignment='left', x=0.01)
            plt.ylabel("norm of absolute position integration error [m]")
            plt.legend(["x", "y", "z"], loc='upper left')
            plt.yscale("log")
            plt.ylim(y_min_pos, y_max_pos)


            plt.subplot(2,1,2)
            plt.plot(t - t[0], allerrors_vel_norm, linewidth=0.75)
            plt.xlabel("t [h]", horizontalalignment='left', x=0.01)
            plt.ylabel("norm of absolute velocity integration error [m]")
            plt.legend(["x", "y", "z"], loc='upper left')
            plt.yscale("log")
            plt.ylim(y_min_vel, y_max_vel)

            plt.tight_layout()
            plt.savefig(dir_plots + body + '_state_history_integration_error'+output_string+'.png')