"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""

# Vehicles
vehicles = ["MESSENGER",
            "BepiColombo"]

# User inputs
no_bodies = 1
no_arcs_v = [11,12]
bodies     = ["Vehicle"]
parameters = []
dependent_variables = []


for i in range(len(vehicles)):

    v = vehicles[i]
    no_arcs = no_arcs_v[i]

    print("for vehicle:", v)

    dir_cpp_output = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/' + v + '/'
    dir_plots = '/home/rens/Documents/PostProcessing_plots/MercuryStateOrbiter/' + v + '/'


    # Plot integration error
    print( " making plots of the integration errors after backward propagation..." )
    import IntegrationError
    IntegrationError.f(dir_cpp_output, dir_plots, bodies, no_arcs)

    # Make correlation heat map
    print(" making heat map of parameter correlations...")
    import HeatMap
    HeatMap.f(dir_cpp_output, dir_plots, parameters, no_arcs)

    # Plot propagated bodies
    print(" making plots of propagated bodies...")
    import PropagatedBodies
    PropagatedBodies.f(dir_cpp_output, dir_plots, bodies, no_arcs)

    # Plot residuals over time
    print(" making plot of the observation residuals and propagated errors...")
    import Residuals
    Residuals.f(dir_cpp_output, dir_plots, v, no_arcs, False)
    Residuals.f(dir_cpp_output, dir_plots, v, no_arcs, True)

    # Plot residuals over time
    print(" making plot of the propagated errors vs true anomaly...")
    import TrueAnomalyVsError
    TrueAnomalyVsError.f(dir_cpp_output, dir_plots)


    # # Plot parameter history
    # print(" making plots of parameter estimation history...")
    # import ParameterHistory
    # ParameterHistory.f(dir_cpp_output, dir_plots, parameters,  bodies)
    #
    # # Plot dependent variable history
    # print(" making plots of dependent variable history...")
    # import DependentVariableHistory
    # DependentVariableHistory.f(dir_cpp_output, dir_plots, dependent_variables)
    #
    # # Make correlation heat map
    # print(" making heat map of parameter correlations...")
    # import HeatMap
    # HeatMap.f(dir_cpp_output, dir_plots, parameters)


print("done!")