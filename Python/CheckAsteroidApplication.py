"""
purpose: checks whether the propagation in "asteroidConsiderCovariance.cpp" is similar to the main application

made on august 12th 2020
"""

def f(dir_output, dir_asteroids, dir_plots):

    import numpy as np
    import matplotlib.pyplot as plt

    # get state history
    filename = "StatePropagationHistoryMercury.dat"
    statehistory_main = np.genfromtxt(dir_output    + filename, delimiter=',')
    statehistory_ast  = np.genfromtxt(dir_asteroids + filename, delimiter=',')

    t_m = statehistory_main[:,0]
    r_m = statehistory_main[:,1:4]
    v_m = statehistory_main[:,4:7]

    t_a = statehistory_ast[:,0]
    r_a = statehistory_ast[:,1:4]
    v_a = statehistory_ast[:,4:7]


    if not np.array_equal(t_m,t_a):
        print(" ERROR, time arrays not equal! moving on")
        print(len(t_m), len(t_a))
        print(t_m)
        print(t_a)
    else:

        error_r = r_a-r_m
        error_r_norm = np.linalg.norm(error_r, axis=1)

        error_v = v_a-v_m
        error_v_norm = np.linalg.norm(error_v, axis=1)

        # print maximum errors
        max_pos_error = np.max(error_r_norm)
        index_pos_error = np.argmax(error_r_norm)
        r_m_norm = np.linalg.norm(r_m, axis=1)
        rel_pos_error = (error_r_norm[index_pos_error] / r_m_norm[index_pos_error])
        print(" max pos error: ", max_pos_error)
        print(" relative: ", rel_pos_error)

        max_vel_error = np.max(error_v_norm)
        index_vel_error = np.argmax(error_v_norm)
        v_m_norm = np.linalg.norm(v_m, axis=1)
        rel_vel_error = (error_v_norm[index_vel_error] / v_m_norm[index_vel_error])
        print(" max pos error: ", max_vel_error)
        print(" relative: ", rel_vel_error)


        #plot errors of all 6 state parameters
        fig1 = plt.figure(figsize=(16, 10))

        plt.subplot(2,1,1)
        for i in range(3):
            plt.plot(t_m, abs(error_r[:,i]))
        plt.xlabel("t [s]")
        plt.ylabel("position difference [m]")
        plt.legend(["x", "y", "z"], loc='upper left')
        plt.yscale("log")

        plt.subplot(2,1,2)
        for i in range(3):
            plt.plot(t_m, abs(error_v[:,i]))
        plt.xlabel("t [s]")
        plt.ylabel("velocity difference [m/s]")
        plt.legend(["x", "y", "z"], loc='upper left')
        plt.yscale("log")

        plt.tight_layout()
        plt.savefig(dir_plots + "asteroidapplicationcheck.png")


        # plot errors of norms
        fig2 = plt.figure(figsize=(16, 10))

        plt.subplot(2, 1, 1)
        plt.plot(t_m, error_r_norm)
        plt.xlabel("t [s]")
        plt.ylabel("position norm difference [m]")
        plt.yscale("log")

        plt.subplot(2, 1, 2)
        plt.plot(t_m, error_v_norm)
        plt.xlabel("t [s]")
        plt.ylabel("velocity norm difference [m/s]")
        plt.yscale("log")

        plt.tight_layout()
        plt.savefig(dir_plots + "asteroidapplicationcheck_norm.png")



