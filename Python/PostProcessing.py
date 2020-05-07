"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""


################
#### INPUTS ####
################

publication_string = ["Genova2018","Schettino2015"]

# Directories
dir_application = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/'

# Bodies included
no_bodies = 1
no_arcs = 1
bodies     = ["Mercury"]

dependent_variables = ["Venus_CG", "Earth_CG", "Moon_CG", "Mars_CG", "Jupiter_CG", "Saturn_CG",
                       "Sun_CG",
                       "exclude", #J1
                       "Sun_J2",]

for ps in publication_string:
    dir_cpp_output = dir_application + 'Output/' + 'Output' + ps + "/"
    dir_plots = '/home/rens/Documents/PostProcessing_plots/thesis_v1/' + ps + "/"
    json_file = dir_application + 'Input/' + 'inputs_' + ps + '.json'

    print(" ")
    print(">>>> FOR INPUTS OF PUBLICATION:", ps)

    # from input file, get which additional things were estimated and add to lists above
    parameters = ["X_Mer", "Y_Mer", "Z_Mer",
                  "Vx_Mer", "Vy_Mer", "Vz_Mer"]

    import json
    print(" json file used as input:" + json_file)
    with open(json_file) as f:
        json_input = json.load(f)

    if json_input["calculateSchwarzschildCorrection"]:
        parameters.append("gamma")
        parameters.append("beta")
        dependent_variables.append("Sun_SS")

    if json_input["includeSEPViolationAcceleration"]:
        if not json_input["useNordtvedtConstraint"]:
            parameters.append("Nordtvedt")
        dependent_variables.append("Sun_SEP")

    if json_input["includeTVGPAcceleration"]:
        parameters.append("TVGP")
        dependent_variables.append("Sun_TVGP")

    parameters.append("mu_Sun")
    parameters.append("J2_Sun")

    print(" ", parameters)

    vehicle = json_input["vehicle"]

    #################
    #### OUTPUTS ####
    #################

    # # Plot integration error
    # print("---- making plots of the integration errors after backward propagation ----" )
    # import IntegrationError
    # IntegrationError.f(dir_cpp_output, dir_plots, bodies, no_arcs)

    # Plot propagated bodies
    print("---- making plots of propagated bodies ----")
    import PropagatedBodies
    PropagatedBodies.f(dir_cpp_output, dir_plots, bodies[0], no_arcs)

    # Plot parameter history
    print("---- making plots of parameter estimation history ----")
    import ParameterHistory
    ParameterHistory.f(dir_cpp_output, dir_plots, parameters, bodies, json_input)

    # Plot residuals over time
    print("---- making plot of the observation residuals and propagated errors ----")
    import Residuals
    Residuals.f(dir_cpp_output, dir_plots, bodies[0], no_arcs, False)

    print("---- making plot of the observation residuals and propagated errors in the RSW frame ----")
    Residuals.f(dir_cpp_output, dir_plots, bodies[0], no_arcs, True)

    print("---- plotting errors interpolated from results of the mercury orbiter ----")
    import InterpolatedErrors
    InterpolatedErrors.f(dir_cpp_output, dir_plots, bodies[0], no_arcs, vehicle)

    # Plot dependent variable history
    print("---- making plots of dependent variable history ----")
    import DependentVariableHistory
    DependentVariableHistory.f(dir_cpp_output, dir_plots, dependent_variables)

    # Make correlation heat map
    print("---- making heat map of parameter correlations ----")
    import HeatMap
    HeatMap.f(dir_cpp_output, dir_plots, parameters, no_arcs)

