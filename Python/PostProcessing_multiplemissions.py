"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""

from os import scandir, mkdir, path
from matplotlib.pyplot import close

################
#### INPUTS ####
################

# Directories output simulations
dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'
dir_output = dir_application + 'Output/'
subdirs = [f.path for f in scandir(dir_output) if f.is_dir()]

# subdirs = [dir_output + "Fienga2019_reality1_estimation1_testCeres"]
print(subdirs)

# Bodies included
bodies = ["Mercury"]
vehicle = "vehicle"
no_bodies = len(bodies)
no_arcs = 1


for s in subdirs:

    if 'reality' not in s:
        continue

    print(" ")
    print(">>>> FOR INPUTS OF PUBLICATION:", s)

    if s[len(dir_output):] != "Xu2017onlyBepi_reality3_estimation1":
        continue

    if s[len(dir_output):] == "Fienga2019_reality1_estimation1_testCeres":
        ps = "Fienga2019"
        reality = 1
        estimation = 1
        continue
    elif s[len(dir_output):] == "Fienga2019_reality1_estimation1_testWithoutFlybys":
        ps = "Fienga2019"
        reality = 1
        estimation = 1
        continue
    elif s[len(dir_output):] == "Fienga2019_reality1_estimation1_onlyUseFirst4Asteroids":
        ps = "Fienga2019"
        reality = 1
        estimation = 1
        continue
    elif s[len(dir_output):] == "VeryLargeAmplitude_reality3_estimation3_testWithoutEstimatingJ2":
        ps = "VeryLargeAmplitude"
        reality = 3
        estimation = 3
    elif s[len(dir_output):] == "Xu2017_reality3_estimation3_testWithoutEstimatingJ2":
        ps = "Xu2017"
        reality = 3
        estimation = 3
    else:
        ps = s[len(dir_output):-21]
        reality = int(s[-13])
        estimation = int(s[-1])


    s += '/'
    json_file = dir_application + 'Input/' + 'inputs_multiplemissions_' + ps + '.json'
    dir_plots = '/home/rens/Documents/PostProcessing_plots/thesis_multiplemissions/' + s[len(dir_output):]
    if not path.exists(dir_plots):
        mkdir(dir_plots)


    #d dependent variables
    dependent_variables = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon",
                           "Sun_CG",
                           "exclude", #J1
                           "Sun_J2"]

    if reality%2 == 0:
        dependent_variables.append("exclude") # J3
        dependent_variables.append("J4_Sun")


    # body state parameters
    parameters = ["X_Mer", "Y_Mer", "Z_Mer",
                  "Vx_Mer", "Vy_Mer", "Vz_Mer"]

    if no_bodies > 1:
        b = 1
        p = parameters
        while b < no_bodies:
            parameters = parameters + p
            b += 1


    #additional parameters
    import json
    print(" json file used as input:" + json_file)
    with open(json_file) as f:
        json_input = json.load(f)

    if s[len(dir_output):] == "Fienga2019_reality1_estimation1_testCeres/":
        parameters.append("mu_Ceres")
        dependent_variables = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon",
                               "Ceres",
                               "Sun_CG",
                               "exclude",  # J1
                               "Sun_J2"]

    if json_input["calculateSchwarzschildCorrection"]:
        if not json_input["gammaIsAConsiderParameter"]:
            parameters.append("gamma")
        parameters.append("beta")
        dependent_variables.append("Sun_SS")
        dependent_variables.append("Sun_SS_α")
        if json_input["calculateLenseThirringCorrection"]:
            dependent_variables.append("Sun_LT")

    if json_input["includeSEPViolationAcceleration"]:
        dependent_variables.append("Sun_SEP")
        if not json_input["useNordtvedtConstraint"]: #comment this line if nordtvedt is enforced in the estimation
            parameters.append("Nordtvedt")

    if json_input["estimatePPNalphas"]:
        parameters.append("alpha1")
        parameters.append("alpha2")

    if json_input["calculateLenseThirringCorrection"] and json_input["estimateSunAngularMomentum"]:
        parameters.append("S_Sun")

    if json_input["includeTVGPAcceleration"]:
        parameters.append("TVGP")
        dependent_variables.append("Sun_TVGP")

    if ps == "AttemptWithMu":
        parameters.append("mu_Sun")

    if estimation > 2:
        parameters.append("J2_A")
        if estimation%2 == 0:
            parameters.append("J4_A")

    parameters.append("J2_Sun")

    if reality % 2 == 0:
        parameters.append("J4_Sun")

    if s[len(dir_output):] == "VeryLargeAmplitude_reality3_estimation3_testWithoutEstimatingJ2/" or \
            s[len(dir_output):] == "Xu2017_reality3_estimation3_testWithoutEstimatingJ2/":
        parameters = ["X_Mer", "Y_Mer", "Z_Mer", "Vx_Mer", "Vy_Mer", "Vz_Mer", "J2_A"]

    print(" ", parameters)


    #################
    #### OUTPUTS ####
    #################

    # Analyse consider covariance due to asteroids
    print("---- analysing consider covariance due to asteroids ----")
    import addedFormalErrorDueToAsteroids
    addedFormalErrorDueToAsteroids.f(s, dir_plots, dir_application, parameters)

    # Make correlation heat map
    print("---- making heat map of parameter correlations ----")
    import HeatMap
    HeatMap.f(s, dir_plots, parameters, no_arcs)

    # Check observation weights
    print("---- checking observation weights ----")
    import CheckWeights
    CheckWeights.f(s, dir_plots)

    # Check consistency asteroid application and main application
    from os.path import isdir
    print(dir_application + 'Input/asteroids_multiplemissions/')
    if isdir(dir_application + 'Input/asteroids_multiplemissions/'):
        print("---- comparing orbits asteroid application and main application----")
        import CheckAsteroidApplication
        CheckAsteroidApplication.f(s, dir_application + 'Input/asteroids_multiplemissions/', dir_plots)

        if s[len(dir_output):] == "Fienga2019_reality1_estimation1_testCeres":
            print("---- making plots of dependent variables of asteroids ----")
            import DependentVariableHistory_Asteroids
            DependentVariableHistory_Asteroids.f(dir_application + 'Input/asteroids_multiplemissions/', dir_plots)

    # Plot dependent variable history
    # print("---- making plots of dependent variable history ----")
    # import DependentVariableHistory
    # DependentVariableHistory.f(s, dir_plots, dependent_variables)

    print("---- making plots of dependent variable history, planets grouped ----")
    import DependentVariableHistory_Grouped
    DependentVariableHistory_Grouped.f(s, dir_plots, dependent_variables)

    # Plot integration error
    print("---- making plots of the integration errors wrt SPICE and backwards integration ----" )
    import IntegrationErrorWrtSPICE
    IntegrationErrorWrtSPICE.f(s, dir_plots, bodies[0], no_arcs)

    close('all')

    # Plot propagated bodies
    print("---- making plots of propagated bodies ----")
    import PropagatedBodies
    PropagatedBodies.f(s, dir_plots, bodies[0], no_arcs)

    # Plot parameter history
    print("---- making plots of parameter estimation history ----")
    import ParameterHistory
    ParameterHistory.f(s, dir_plots, parameters, no_bodies, json_input, False)

    # Plot residuals over time
    print("---- making plot of the observation residuals and propagated errors ----")
    import Residuals
    Residuals.f(s, dir_plots, bodies[0], no_arcs, False)

    print("---- making plot of the observation residuals and propagated errors in the RSW frame ----")
    Residuals.f(s, dir_plots, bodies[0], no_arcs, True)

    print("---- plotting errors interpolated from results of the mercury orbiter ----")
    import InterpolatedErrors
    InterpolatedErrors.f(s, dir_plots, bodies[0], no_arcs, vehicle)

    close('all')

