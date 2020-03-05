"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""


# Directories
dir_cpp_output = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

# User inputs
no_bodies = 1
bodies     = ["Mercury"]
parameters = ["X_Mer", "Y_Mer", "Z_Mer",
              "Vx_Mer", "Vy_Mer", "Vz_Mer",
              "PPN_beta",
              "mu_Sun",
              "J2_Sun",]


# Plot propagated bodies
print("making plots of propagated bodies...")
import PropagatedBodies
PropagatedBodies.f(dir_cpp_output, dir_plots, parameters, bodies)

# Plot parameter history
print("making plots of parameter estimation history...")
import ParameterHistory
ParameterHistory.f(dir_cpp_output, dir_plots, parameters,  bodies)

# Make correlation heat map
print("making heat map of parameter correlations...")
import HeatMap
HeatMap.f(dir_cpp_output, dir_plots, parameters)


print("done!")