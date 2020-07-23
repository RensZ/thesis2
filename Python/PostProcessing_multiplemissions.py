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

    if json_input["calculateSchwarzschildCorrection"]:
        if not json_input["gammaIsAConsiderParameter"]:
            parameters.append("gamma")
        parameters.append("beta")
        dependent_variables.append("Sun_SS")
        dependent_variables.append("Sun_SS_α")
        if json_input["calculateLenseThirringCorrection"]:
            dependent_variables.append("Sun_LT")

    if json_input["includeSEPViolationAcceleration"]:
        # if not json_input["useNordtvedtConstraint"]:
        parameters.append("Nordtvedt")
        dependent_variables.append("Sun_SEP")

    if json_input["estimatePPNalphas"]:
        parameters.append("alpha1")
        parameters.append("alpha2")

    if json_input["calculateLenseThirringCorrection"] and json_input["estimateSunAngularMomentum"]:
        parameters.append("S_Sun")

    if json_input["includeTVGPAcceleration"]:
        parameters.append("TVGP")
        dependent_variables.append("Sun_TVGP")

    parameters.append("mu_Sun")

    if estimation > 2:
        parameters.append("J2_A")
        if estimation%2 == 0:
            parameters.append("J4_A")

    parameters.append("J2_Sun")
    if reality % 2 == 0:
        parameters.append("J4_Sun")

    print(" ", parameters)


    #################
    #### OUTPUTS ####
    #################

    # Plot dependent variable history
    print("---- making plots of dependent variable history ----")
    import DependentVariableHistory
    DependentVariableHistory.f(s, dir_plots, dependent_variables)

    # Plot integration error
    print("---- making plots of the integration errors wrt SPICE and backwards integration ----" )
    import IntegrationErrorWrtSPICE
    IntegrationErrorWrtSPICE.f(s, dir_plots, bodies[0], no_arcs)

    # Plot propagated bodies
    print("---- making plots of propagated bodies ----")
    import PropagatedBodies
    PropagatedBodies.f(s, dir_plots, bodies[0], no_arcs)

    # Plot parameter history
    print("---- making plots of parameter estimation history ----")
    import ParameterHistory
    ParameterHistory.f(s, dir_plots, parameters, no_bodies, json_input, False)

    # Make correlation heat map
    print("---- making heat map of parameter correlations ----")
    import HeatMap
    HeatMap.f(s, dir_plots, parameters, no_arcs)

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

