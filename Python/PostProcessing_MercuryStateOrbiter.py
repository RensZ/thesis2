"""
Created by Rens van der Zwaard on 2020-2-25

Purpose: wrapper file for all the post-processing of thesis_v1.cpp

"""


# Directories
dir_cpp_output = '/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/mercuryOrbiterState/MercuryOrbiterStateEstimation/'
dir_plots = '/home/rens/Documents/PostProcessing_plots/'

# User inputs
no_bodies = 1
bodies     = ["Vehicle"]
parameters = ["X_Veh", "Y_Veh", "Z_Veh",
              "Vx_Veh", "Vy_Veh", "Vz_Veh",
              "PPN_gamma",
              "CR_Veh",
              "J2_Sun"]
dependent_variables = ["Venus_CG",
                       "Earth_CG",
                       "Moon_CG",
                       "Mars_CG",
                       "Jupiter_CG",
                       "Saturn_CG",
                       "Sun_CG",
                       "exclude",
                       "Sun_J2",
                       "Sun_SS",
                       "Sun_TVGP"]


# Plot propagated bodies
print("making plots of propagated bodies...")
import PropagatedBodies
PropagatedBodies.f(dir_cpp_output, dir_plots, parameters, bodies)

# Plot residuals over time
print("making plot of the observation residuals...")
import Residuals
Residuals.f(dir_cpp_output, dir_plots, bodies)

# # Plot parameter history
# print("making plots of parameter estimation history...")
# import ParameterHistory
# ParameterHistory.f(dir_cpp_output, dir_plots, parameters,  bodies)
#
# # Plot dependent variable history
# print("making plots of dependent variable history...")
# import DependentVariableHistory
# DependentVariableHistory.f(dir_cpp_output, dir_plots, dependent_variables)
#
# # Make correlation heat map
# print("making heat map of parameter correlations...")
# import HeatMap
# HeatMap.f(dir_cpp_output, dir_plots, parameters)


print("done!")