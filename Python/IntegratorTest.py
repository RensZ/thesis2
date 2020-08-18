"""
purpose: to process the output of integratorTest.cpp, to compare integrator settings

"""

from os import scandir
import matplotlib.pyplot as plt
import numpy as np

integrators = ["ABM","RK4","RK7"]

plotsteps = [300, 900, 1800, 3600, 7200, 14400]

dir_plots = '/home/rens/Documents/PostProcessing_plots/integratorTest/'
dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_output = dir_application + 'Output/integratorTest/'
allfiles = [f.path for f in scandir(dir_output) ]

entries1_timesteps = []
entries1_errors_vs_timestep = []
legend_entries1_errors_vs_timestep = []
entries2_timesteps = []
entries2_errors_vs_timestep = []
legend_entries2_errors_vs_timestep = []

for integratorString in np.asarray(integrators):

    timesteps = []
    timearrays = []
    posarrays = []
    backarrays = []
    max_errors = []
    max_backerrors = []

    for f in allfiles:

        #only continue if its a state propagation history file of the current integrator
        if len(f) > 120 or integratorString not in f:
            continue

        #find the timestep
        if len(f) == 118:
            timestep = f[-7:-4]
        elif len(f) == 119:
            timestep = f[-8:-4]
        elif len(f) == 120:
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


    #sort arrays (makes plots look better)
    sorted_indices = np.argsort(timesteps)
    timesteps = np.asarray(timesteps)[sorted_indices]
    timearrays = np.asarray(timearrays)[sorted_indices]
    posarrays = np.asarray(posarrays)[sorted_indices]
    backarrays = np.asarray(backarrays)[sorted_indices]
    max_backerrors = np.asarray(max_backerrors)[sorted_indices]


    #plot differences forward and backward propagation
    print("plotting forward - backward propagation")
    legend_entries = []
    fig1 = plt.figure(figsize=(12,6))

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
    plt.savefig(dir_plots + "forward_minus_backwards_" + integratorString + ".png")


    #plot differences, assuming the integration with minimum timestep is the truth
    print("plotting state history - benchmark")
    min_timestep = int(np.argmin(timesteps))
    min_timearray = timearrays[min_timestep]
    min_posarray = posarrays[min_timestep]

    legend_entries = []
    fig1 = plt.figure(figsize=(12,6))

    for i in range(len(timesteps)):

        if i == min_timestep:
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

        max_errors.append(max(norm_difference))

        if timesteps[i] in plotsteps:

            plt.plot(current_timearray_matched, norm_difference, linewidth=0.75)
            legend_entries.append(str(timesteps[i]))

    plt.xlabel("time [s]")
    plt.ylabel("error norm [m]")
    plt.yscale("log")
    plt.legend(legend_entries)
    plt.tight_layout()
    plt.savefig(dir_plots + "benchmark_error_" + integratorString + ".png")

    # entries for the plot after this loop
    entries1_timesteps.append(timesteps[1:])
    entries1_errors_vs_timestep.append(max_errors)
    legend_entries1_errors_vs_timestep.append(integratorString)

    entries2_timesteps.append(timesteps)
    entries2_errors_vs_timestep.append(max_backerrors)
    legend_entries2_errors_vs_timestep.append(integratorString)


# plotting error vs timestep
print("plotting max errors vs timestep")

fig2 = plt.figure(figsize=(12,8))

plt.subplot(1,2,1)
for e in range(len(entries1_errors_vs_timestep)):
    plt.plot(entries1_timesteps[e], entries1_errors_vs_timestep[e], linewidth=0.5, marker='o', markersize=5)

plt.xlabel('time step [s]')
plt.ylabel('maximum error norm w.r.t. benchmark [m]')
plt.xlim(np.min(timesteps)-10.0, np.max(timesteps)+1000.0)
plt.xscale('log')
plt.yscale('log')
plt.legend(legend_entries1_errors_vs_timestep)

plt.subplot(1,2,2)
for e in range(len(entries2_errors_vs_timestep)):
    plt.plot(entries2_timesteps[e], entries2_errors_vs_timestep[e], linewidth=0.5, marker='o', markersize=5)

plt.xlabel('time step [s]')
plt.ylabel('maximum error forward-backward [m]')
plt.xlim(np.min(timesteps)-10.0, np.max(timesteps)+1000.0)
plt.xscale('log')
plt.yscale('log')
plt.legend(legend_entries2_errors_vs_timestep)


plt.tight_layout()
plt.savefig(dir_plots + "maxerror_vs_timestep.png")


plt.close('all')

print(np.max(posarrays[0]))
print(np.max(posarrays[0])*1E-16)