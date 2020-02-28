"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""


# Inputs
dir_cpp_output = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

no_bodies = 1
parameters = ["X_Mer", "Y_Mer", "Z_Mer",
              "Vx_Mer", "Vy_Mer", "Vz_Mer",
              "PPN_gamma",
              "PPN_beta",
              "mu_Sun",
              "J2_Sun"]


# Plot parameter history
print("making plot of parameter history...")
import ParameterHistory
ParameterHistory.f(dir_cpp_output, parameters, dir_plots, no_bodies)

# Make correlation heat map
print("making heat map of parameter correlations...")
import HeatMap
HeatMap.f(dir_cpp_output, parameters, dir_plots)

