"""
purpose: to process the output of integratorTest.cpp, to compare integrator settings

"""

from os import scandir
import matplotlib.pyplot as plt
import numpy as np

integrators = ["RK4","RK7","ABM","A12"]
plotsteps = [300, 600, 900, 1800, 3600, 7200, 14400]
# outputSubFolders = ["integratorTest/", "integratorTest_long/"]
outputSubFolders = ["integratorTest_long/"]
dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'

for o in outputSubFolders:

    dir_plots = '/home/rens/Documents/PostProcessing_plots/' + o
    dir_output = dir_application + 'Output/' + o
    allfiles = [f.path for f in scandir(dir_output) ]
    print(allfiles)

    entries1_timesteps = []
    entries1_lengths = []
    entries1_errors_vs_timestep = []
    legend_entries1_errors_vs_timestep = []
    entries2_timesteps = []
    entries2_lengths = []
    entries2_errors_vs_timestep = []
    legend_entries2_errors_vs_timestep = []

    # get benchmark data
    minString = "RK4"
    if o == "integratorTest":
        min_timestep = 750
    else:
        min_timestep = 600
    min_data = np.genfromtxt(dir_output + "StatePropagationHistory"+minString + str(min_timestep) + ".dat", delimiter=',')
    # min_back = np.genfromtxt(dir_output + "StatePropagationHistoryBackwardsRK436000.dat", delimiter='')
    min_timearray = min_data[:,0]
    min_posarray = min_data[:,1:4]

    for integratorString in np.asarray(integrators):

        extralength = 0
        if o == "integratorTest_long/":
            extralength = 5

        if integratorString == "A12":
            saveString = "ABM12"
        elif integratorString == "ABM":
            saveString = "ABM8"
        elif integratorString == "RK7":
            saveString = "RK78"
        else:
            saveString = integratorString

        timesteps = []
        timearrays = []
        posarrays = []
        backarrays = []
        lengtharrays = []
        max_errors = []
        max_backerrors = []

        for f in allfiles:

            #only continue if its a state propagation history file of the current integrator
            if len(f) > 120+extralength or integratorString not in f:
                continue

            #find the timestep
            if len(f) == 118+extralength:
                timestep = f[-7:-4]
            elif len(f) == 119+extralength:
                timestep = f[-8:-4]
            elif len(f) == 120+extralength:
                timestep = f[-9:-4]
            else:
                timestep = 0
                print("ERROR: file length not supported")

            print("getting data for integrator " + integratorString + " and timestep", timestep)
            timesteps.append(float(timestep))

            # get names all 3 files
            f_b = dir_output + "StatePropagationHistoryBackwards" + integratorString + timestep + ".dat"
            data_b = np.genfromtxt(f_b, delimiter=',')
            # f_s = dir_output + "spiceStatesAtPropagationTimesABM" + timestep + ".dat"

            # get data
            statehistory = np.genfromtxt(f, delimiter=',')
            timeforward = statehistory[:,0]
            timebackward = data_b[:,0]
            if not np.array_equal(timeforward, timebackward):
                print("ERROR: time arrays not equal!")
                continue

            timearrays.append( timeforward )
            posarrays.append( statehistory[:,1:4] )
            backarrays.append( posarrays[-1] - data_b[:,1:4] )
            max_backerrors.append( np.max(np.linalg.norm(backarrays[-1], axis=1)))
            lengtharrays.append( len(timeforward) )


        #sort arrays (makes plots look better)
        sorted_indices = np.argsort(timesteps)
        timesteps = np.asarray(timesteps)[sorted_indices]
        timearrays = np.asarray(timearrays)[sorted_indices]
        posarrays = np.asarray(posarrays)[sorted_indices]
        backarrays = np.asarray(backarrays)[sorted_indices]
        max_backerrors = np.asarray(max_backerrors)[sorted_indices]
        lengtharrays = np.asarray(lengtharrays)[sorted_indices]


        #plot differences forward and backward propagation
        print("plotting forward - backward propagation")
        legend_entries = []
        fig1 = plt.figure(figsize=(10,5))

        for i in range(len(timesteps)):
            if timesteps[i] in plotsteps:

                current_timearray = timearrays[i]
                current_backarray = backarrays[i]

                plt.plot(current_timearray, np.linalg.norm(current_backarray, axis=1), linewidth=0.75)
                legend_entries.append(str(timesteps[i]))

        plt.xlabel("time [s]")
        plt.ylabel("error norm [m]")
        plt.yscale("log")
        plt.legend(legend_entries)
        plt.tight_layout()
        plt.savefig(dir_plots + "forward_minus_backwards_" + saveString + ".png")


        #plot differences, assuming the integration with minimum timestep is the truth
        print("plotting state history - benchmark")

        legend_entries = []
        fig1 = plt.figure(figsize=(10,5))

        for i in range(len(timesteps)):

            if timesteps[i] == min_timestep and integratorString == minString:
                continue

            current_timearray = timearrays[i]
            current_posarray = posarrays[i]
            current_backarray = backarrays[i]

            matching_indices_1 = np.where(np.in1d(min_timearray, current_timearray))[0]
            min_timearray_matched = min_timearray[matching_indices_1]
            min_posarray_matched = min_posarray[matching_indices_1,:]

            matching_indices_2 = np.where(np.in1d(current_timearray, min_timearray_matched))[0]
            current_timearray_matched = current_timearray[matching_indices_2]
            current_posarray_matched = current_posarray[matching_indices_2]

            pos_difference = current_posarray_matched - min_posarray_matched
            norm_difference = np.linalg.norm(pos_difference, axis=1)

            #print(timesteps[i], max(norm_difference))
            max_errors.append(max(norm_difference))

            if timesteps[i] in plotsteps:

                plt.plot(current_timearray_matched, norm_difference, linewidth=0.75)
                legend_entries.append(str(timesteps[i]))


        plt.xlabel("time [s]")
        plt.ylabel("error norm [m]")
        plt.yscale("log")
        plt.legend(legend_entries)
        plt.tight_layout()
        plt.savefig(dir_plots + "benchmark_error_" + saveString + ".png")

        # entries for the plot after this loop
        if integratorString == minString:
            entries1_timesteps.append(timesteps[np.where(timesteps != min_timestep)[0]])
            entries1_lengths.append(lengtharrays[np.where(timesteps != min_timestep)[0]])
        else:
            entries1_timesteps.append(timesteps)
            entries1_lengths.append(lengtharrays)
        entries1_errors_vs_timestep.append(max_errors)
        legend_entries1_errors_vs_timestep.append(saveString)

        entries2_timesteps.append(timesteps)
        entries2_lengths.append(lengtharrays)
        entries2_errors_vs_timestep.append(max_backerrors)
        legend_entries2_errors_vs_timestep.append(saveString)

    # where is the forward - backward error minimum?
    print("min error of each integrator:")
    for e in range(len(entries2_errors_vs_timestep)):
        alltimes = entries1_timesteps[e]
        allerrors = entries1_errors_vs_timestep[e]
        print(allerrors)
        print(integrators[e], alltimes[np.argmin(allerrors)], min(allerrors))

    # plotting error vs timestep
    print("plotting max errors vs timestep")

    fig2 = plt.figure(figsize=(10,5))

    plt.subplot(1,2,1)
    for e in range(len(entries1_errors_vs_timestep)):
        plt.plot(entries1_timesteps[e], entries1_errors_vs_timestep[e], linewidth=0.5, marker='o', markersize=5)

    plt.xlabel('time step [s]')
    plt.ylabel('maximum error norm w.r.t. benchmark [m]')
    plt.xlim(np.min(timesteps)-100.0, np.max(timesteps)+10000.0)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(legend_entries1_errors_vs_timestep)

    plt.subplot(1,2,2)
    for e in range(len(entries2_errors_vs_timestep)):
        plt.plot(entries2_timesteps[e], entries2_errors_vs_timestep[e], linewidth=0.5, marker='o', markersize=5)

    plt.xlabel('time step [s]')
    plt.ylabel('maximum error forward-backward [m]')
    plt.xlim(np.min(timesteps)-100.0, np.max(timesteps)+10000.0)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(legend_entries2_errors_vs_timestep)


    plt.tight_layout()
    plt.savefig(dir_plots + "maxerror_vs_timestep.png")




    # plotting error vs timestep
    print("plotting max errors vs function evaluations")

    fig3 = plt.figure(figsize=(10,5))

    plt.subplot(1,2,1)
    for e in range(len(entries1_errors_vs_timestep)):
        if legend_entries1_errors_vs_timestep[e] == "RK78":
            funceval = 13
        else:
            funceval = 4
        numberOfEvaluations = entries1_lengths[e] * funceval
        plt.plot(numberOfEvaluations, entries1_errors_vs_timestep[e], linewidth=0.5, marker='o', markersize=5)

        print("From backward integration, the 10 options with the lowest error for integrator " + legend_entries1_errors_vs_timestep[e])
        sorted_indices1 = np.argsort(entries2_errors_vs_timestep[e])[0:10]
        t_sorted = np.asarray(entries2_timesteps[e])[sorted_indices1]
        n_sorted = np.asarray(numberOfEvaluations)[sorted_indices1]
        e_sorted = np.asarray(entries2_errors_vs_timestep[e])[sorted_indices1]
        print(t_sorted)
        print(n_sorted)
        print(e_sorted)


    plt.xlabel('function evaluations')
    plt.ylabel('maximum error norm w.r.t. benchmark [m]')
    plt.xlim(np.min(numberOfEvaluations)*0.95, np.max(numberOfEvaluations)+1.05)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(legend_entries1_errors_vs_timestep)

    print(" ---- ")
    plt.subplot(1,2,2)
    for e in range(len(entries2_errors_vs_timestep)):
        if legend_entries2_errors_vs_timestep[e] == "RK78":
            funceval = 13
        else:
            funceval = 4
        numberOfEvaluations = entries2_lengths[e]*funceval
        plt.plot(numberOfEvaluations, entries2_errors_vs_timestep[e], linewidth=0.5, marker='o', markersize=5)

        print("From benchmark data, the 10 options with the lowest error for integrator " + legend_entries2_errors_vs_timestep[e])
        sorted_indices2 = np.argsort(entries2_errors_vs_timestep[e])[0:10]
        t_sorted = np.asarray(entries2_timesteps[e])[sorted_indices2]
        n_sorted = np.asarray(numberOfEvaluations)[sorted_indices2]
        e_sorted = np.asarray(entries2_errors_vs_timestep[e])[sorted_indices2]
        print(t_sorted)
        print(n_sorted)
        print(e_sorted)



    plt.xlabel('function evaluations')
    plt.ylabel('maximum error forward-backward [m]')
    plt.xlim(np.min(numberOfEvaluations)*0.95, np.max(numberOfEvaluations)+1.05)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(legend_entries2_errors_vs_timestep)

    plt.tight_layout()
    plt.savefig(dir_plots + "maxerror_vs_functionevaluations.png")

    plt.close('all')