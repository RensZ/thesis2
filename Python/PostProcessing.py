"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""

from os import mkdir, path
from matplotlib.pyplot import close

################
#### INPUTS ####
################

publication_string = [# "Genova2018_testGMSunAsConsiderParameter",
                      "Schettino2015_enforceNordtvedtConstraint",
                      "Imperi2018_nvtrue_enforceNordtvedtConstraint",
                      "Imperi2018_nvtrue_flybys_alphas_enforceNordtvedtConstraint",
                      "Genova2018_enforceNordtvedtConstraint",
                      "Schettino2015",
                      "Imperi2018_nvtrue",
                      "Imperi2018_nvfalse",
                      "Genova2018",
                      "Schettino2015_alphas",
                      "Imperi2018_nvtrue_flybys_alphas",
                      "Imperi2018_nvfalse_flybys_alphas",
]

#"MESSENGER_and_BepiColombo",
#"MESSENGER_and_BepiColombo_timevariableJ2",
# "Schettino2015_testCeres",
# "Genova2018_testCeres",
# "Genova2018_testGamma",
# "Imperi2018_nvtrue_flybys",
# "Imperi2018_nvfalse_flybys",
# "Schettino2015_alphas",

# Directories
dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'

# Bodies included
bodies = ["Mercury"]
no_arcs = 1



for ps in publication_string:

    print(" ")
    print(">>>> FOR INPUTS OF PUBLICATION:", ps)

    if ps == "MESSENGER_and_BepiColombo_multiarc":
        no_bodies = 2
        ps_json = "MESSENGER_and_BepiColombo"
    else:
        no_bodies = 1
        ps_json = ps

    enforceNordtvedtConstraint = False
    if ps == "Schettino2015_enforceNordtvedtConstraint":
        ps_json = "Schettino2015"
        enforceNordtvedtConstraint = True
    if ps == "Imperi2018_nvtrue_enforceNordtvedtConstraint":
        ps_json = "Imperi2018_nvtrue"
        enforceNordtvedtConstraint = True
    if ps == "Imperi2018_nvtrue_flybys_alphas_enforceNordtvedtConstraint":
        ps_json = "Imperi2018_nvtrue_flybys_alphas"
        enforceNordtvedtConstraint = True
    if ps == "Genova2018_enforceNordtvedtConstraint":
        ps_json = "Genova2018"
        enforceNordtvedtConstraint = True

    if ps == "Genova2018_testGamma" or ps == "Genova2018_testCeres" or ps == "Genova2018_testGMSunAsConsiderParameter":
        ps_json = "Genova2018"
    if ps == "Schettino2015_testCeres" or ps == "Schettino2015_testGMSunAsConsiderParameter":
        ps_json = "Schettino2015"

    dir_cpp_output = dir_application + 'Output/' + ps + "/"
    json_file = dir_application + 'Input/' + 'inputs_' + ps_json + '.json'
    dir_plots = '/home/rens/Documents/PostProcessing_plots/thesis_v1/' + ps + "/"
    if not path.exists(dir_plots):
        mkdir(dir_plots)


    # from input file, get which additional things were estimated and add to lists above
    parameters = ["X_Mer", "Y_Mer", "Z_Mer",
                  "Vx_Mer", "Vy_Mer", "Vz_Mer"]

    if no_bodies > 1:
        b = 1
        p = parameters
        while b < no_bodies:
            parameters = parameters + p
            b += 1

    dependent_variables = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon",
                           "Sun_CG",
                           "exclude",  # J1
                           "Sun_J2"]

    import json
    print(" json file used as input:" + json_file)
    with open(json_file) as f:
        json_input = json.load(f)

    # if json_input["estimateJ4Amplitude"] or json_input["estimateJ4Amplitude"] or json_input["estimateJ4Amplitude"]:
    # dependent_variables.append("exclude") #J3
    # dependent_variables.append("J4_Sun")

    if json_input["calculateSchwarzschildCorrection"]:

        if ps == "Genova2018_testGamma":
            parameters.append("gamma")

        if ps == "Genova2018_testCeres" or ps == "Schettino2015_testCeres":
            parameters.append("mu_Ceres")
            dependent_variables = ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon", "Ceres",
                                   "Sun_CG",
                                   "exclude",  # J1
                                   "Sun_J2"]

        if not json_input["gammaIsAConsiderParameter"]:
            parameters.append("gamma")
        parameters.append("beta")
        dependent_variables.append("Sun_SS")
        dependent_variables.append("Sun_SS_α")
        if json_input["calculateLenseThirringCorrection"]:
            dependent_variables.append("Sun_LT")
        # if json_input["calculateDeSitterAcceleration"]:
        #     dependent_variables.append("Sun_DS")

    if json_input["includeSEPViolationAcceleration"]:
        dependent_variables.append("Sun_SEP")
        if json_input["useNordtvedtConstraint"] == False or enforceNordtvedtConstraint == True: #comment this line if nordtvedt is enforced in the estimation
            parameters.append("Nordtvedt")


    if json_input["estimatePPNalphas"]:
        parameters.append("alpha1")
        parameters.append("alpha2")

    if json_input["calculateLenseThirringCorrection"] and json_input["estimateSunAngularMomentum"]:
        parameters.append("S_Sun")

    if json_input["includeTVGPAcceleration"]:
        parameters.append("TVGP")
        dependent_variables.append("Sun_TVGP")

    if ps != "Genova2018_testGMSunAsConsiderParameter" and ps != "Schettino2015_testGMSunAsConsiderParameter":
        parameters.append("mu_Sun")

    if json_input["estimateJ2Amplitude"]:
        parameters.append("J2_A")
    if json_input["estimateJ2Period"]:
        parameters.append("J2_P")
    if json_input["estimateJ2Phase"]:
        parameters.append("J2_phi")
    if json_input["estimateJ4Amplitude"]:
        parameters.append("J4_A")
    if json_input["estimateJ4Period"]:
        parameters.append("J4_P")
    if json_input["estimateJ4Phase"]:
        parameters.append("J4_phi")

    parameters.append("J2_Sun")
    # # if json_input["estimateJ4Amplitude"] or json_input["estimateJ4Amplitude"] or json_input["estimateJ4Amplitude"]:
    # parameters.append("J4_Sun")


    print(" ", parameters)

    if ps == "MESSENGER_and_BepiColombo" or ps == "MESSENGER_and_BepiColombo_multiarc" or ps == "MESSENGER_and_BepiColombo_timevariableJ2":
        vehicle = "vehicle"
    else:
        vehicle = json_input["vehicle"]

    #################
    #### OUTPUTS ####
    #################

    # # Check consistency asteroid application and main application
    # from os.path import isdir
    # if isdir(dir_cpp_output[:-1]+"_asteroids/"):
    #     print("---- comparing orbits asteroid application and main application----")
    #     import CheckAsteroidApplication
    #     CheckAsteroidApplication.f(dir_cpp_output, dir_cpp_output[:-1] + "_asteroids/", dir_plots)
    #
    #     print("---- making plots of dependent variables of asteroids ----")
    #     import DependentVariableHistory_Asteroids
    #     DependentVariableHistory_Asteroids.f(dir_cpp_output[:-1] + "_asteroids/", dir_plots)

    # Plot dependent variable history
    print("---- making plots of dependent variable history ----")
    import DependentVariableHistory
    DependentVariableHistory.f(dir_cpp_output, dir_plots, dependent_variables)

    print("---- making plots of dependent variable history, planets grouped ----")
    import DependentVariableHistory_Grouped
    DependentVariableHistory_Grouped.f(dir_cpp_output, dir_plots, dependent_variables)

    # Check observation weights
    print("---- checking observation weights ----")
    import CheckWeights
    CheckWeights.f(dir_cpp_output, dir_plots)

    close('all')

    # Plot integration error
    print("---- making plots of the integration errors wrt SPICE and backwards integration ----" )
    import IntegrationErrorWrtSPICE
    IntegrationErrorWrtSPICE.f(dir_cpp_output, dir_plots, bodies[0], no_arcs)

    # Plot propagated bodies
    print("---- making plots of propagated bodies ----")
    import PropagatedBodies
    PropagatedBodies.f(dir_cpp_output, dir_plots, bodies[0], no_arcs)

    # Plot parameter history
    print("---- making plots of parameter estimation history ----")
    import ParameterHistory
    ParameterHistory.f(dir_cpp_output, dir_plots, parameters, no_bodies, json_input, True)

    close('all')

    # Make correlation heat map
    print("---- making heat map of parameter correlations ----")
    import HeatMap
    HeatMap.f(dir_cpp_output, dir_plots, parameters, no_arcs)

    # Plot residuals over time
    print("---- making plot of the observation residuals and propagated errors ----")
    import Residuals
    Residuals.f(dir_cpp_output, dir_plots, bodies[0], no_arcs, False)

    print("---- making plot of the observation residuals and propagated errors in the RSW frame ----")
    Residuals.f(dir_cpp_output, dir_plots, bodies[0], no_arcs, True)

    print("---- plotting errors interpolated from results of the mercury orbiter ----")
    import InterpolatedErrors
    InterpolatedErrors.f(dir_cpp_output, dir_plots, bodies[0], no_arcs, vehicle)


    close('all')

