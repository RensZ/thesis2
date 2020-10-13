import numpy as np
import matplotlib.pyplot as plt

dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_output = dir_application + 'Output/AmplitudeAnalysis/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'
publication = "PaperInputs"

amplitudes = np.genfromtxt("/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/AmplitudeInputs.txt")
# amplitudes = [2.5e-7, 1.0e-7,
#               7.5e-8, 5.0e-8, 2.5e-8, 1.0e-8,
#               7.5e-9, 5.0e-9, 2.5e-9, 1.0e-9,
#               7.5e-10, 5.0e-10, 2.5e-10, 1.0e-10]

amplitudes_relative = np.asarray(amplitudes)/2.25e-7

parameters = ["gamma", "beta", "TVGP", "J2_Sun", "eta"]
parameter_labels = [r"$\gamma$", r"$\beta$", r"$\dot{Gm_{\odot}}/Gm_{\odot}[y^{-1}]$", r"$J_{2\odot}$", r"$\eta$"]
missions = ["onlyMess","onlyBepi",""]
missions_legend = ["only MESSENGER", "only BepiColombo", "Combined"]

gammaIndex = 0
betaIndex = 1
alpha1sigma = 1.9E-5
alpha2sigma = 8.0E-10

no_m = len(missions)
no_a = len(amplitudes)
no_p = len(parameters)

trueErrorArray = np.zeros((no_m, no_a, no_p))
formalErrorArray = np.zeros((no_m, no_a, no_p))

for m in range(no_m):

    for a in range(no_a):

        subfolder = dir_output+publication+missions[m]+"_reality3_estimation1_amp"+str(a)+"/"
        truthParameters = np.genfromtxt(subfolder+"TruthParameters.dat")
        trueErrors = np.genfromtxt(subfolder+"ObservationTrueEstimationError.dat")
        formalErrors = np.genfromtxt(subfolder+"ObservationFormalEstimationErrorWithConsiderIncludingAsteroidsParameters.dat")
        covariance = np.genfromtxt(subfolder+"ConsiderIncludingAsteroidsCovarianceMatrix.dat")

        for p in range(no_p):
            if parameters[p] == "eta":
                trueErrorEta = 4.0*(1.0+trueErrors[6+betaIndex]) - (1.0+trueErrors[6+gammaIndex]) - 3.0
                formalCovarianceEta = 16.0*covariance[6+betaIndex,6+betaIndex] \
                                    + covariance[6+gammaIndex,6+gammaIndex] \
                                    + alpha1sigma**2 \
                                    + ((2.0/3.0)**2)*alpha2sigma**2 \
                                    - 2.0*4.0*covariance[6+gammaIndex,6+betaIndex]
                trueErrorArray[m,a,p] = trueErrorEta
                formalErrorArray[m,a,p] = np.sqrt(formalCovarianceEta)
            else:
                trueErrorArray[m,a,p] = trueErrors[6+p]
                formalErrorArray[m,a,p] = formalErrors[6+p]


fig = plt.figure(figsize=(6,10))
for p in range(no_p):

    plt.subplot(5, 1, p+1)

    for m in range(no_m):
        plt.errorbar(amplitudes_relative, np.abs(trueErrorArray[m,:,p]), xerr=None, yerr=formalErrorArray[m,:,p],
                     fmt='--o', linewidth=0.75, markersize=3.0)
        print(trueErrorArray[m,:,p])
        print(formalErrorArray[m, :, p])
    plt.xscale("log")
    plt.yscale("log")
    if p == no_p-1:
        plt.xlabel(r"Amplitude of $J_{2\odot}$ [relative to $J_{2\odot}$]")
    plt.ylabel("true error "+parameter_labels[p])
    if p == 0:
        plt.legend(missions_legend)
    plt.grid()

plt.tight_layout()
plt.savefig(dir_plots+"parameter_errors_vs_amplitude.png")
plt.close()
