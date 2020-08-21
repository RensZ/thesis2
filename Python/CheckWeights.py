
def f(dir_output, dir_plots):

    import numpy as np
    import matplotlib.pyplot as plt

    #get data
    weightDiagonal = np.genfromtxt(dir_output+"ObservationWeightDiagonal.dat")
    observationTime = np.genfromtxt(dir_output+"interpolatedErrorMatrix.dat")[:,0]
    finalResiduals = np.genfromtxt(dir_output+"ResidualHistory.dat")[:,-1]

    #plot weights vs residuals
    fig = plt.figure(figsize=(16,10))
    plt.plot(weightDiagonal, abs(finalResiduals), "ro", markersize=1)
    plt.xscale("log")
    plt.xlabel("weight")
    plt.yscale("log")
    plt.ylabel("absolute residual [m]")
    plt.grid()
    plt.tight_layout()
    plt.savefig(dir_plots + "weights_vs_residuals.png")

    #plot weights and residuals vs time
    fig2 = plt.figure(figsize=(16,10))

    plt.subplot(2,1,1)
    plt.plot(observationTime, abs(finalResiduals), "ro", markersize=1)
    plt.xlabel("time [s]")
    plt.yscale("log")
    plt.ylabel("absolute residuals [m]")
    plt.grid()

    plt.subplot(2,1,2)
    plt.plot(observationTime, weightDiagonal, "bo", markersize=1)
    plt.xlabel("time [s]")
    plt.yscale("log")
    plt.ylabel("weight")
    plt.grid()

    plt.tight_layout()
    plt.savefig(dir_plots + "weights_and_residuals_vs_time.png")