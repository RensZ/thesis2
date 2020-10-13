"""
Purpose: to test what the effect is of asteroid mass uncertainty

"""

# Park et al 2017, equation (3)
def precessionDueToJ2(J2):

    from math import pi, sin
    R_Sun = 695700E3 #m, nasa fact sheet
    a_Mer = 57.91E9 #m, nasa fact sheet
    e_Mer = 0.2056 #-, nasa fact sheet
    i_Mer = 3.4*pi/180.0 #rad, lit study
    P_Mer = 87.97*24.0*3600.0 #s, google
    century = 100.0*365.25*24.0*3600.0 #s

    n_Mer = 2.0*pi / P_Mer #rad/s

    precessionPerSec = 1.5 \
        * ( n_Mer * J2 / ((1.0 - e_Mer**2)**2.0) ) \
        * ( (R_Sun / a_Mer)**2.0 ) \
        * ( 1.0 - 1.5*(sin(i_Mer)**2.0) )

    precessionPerCentury = century*precessionPerSec

    return precessionPerCentury

# Will 2018 eq (1) term 2
def precessionDueToCeres(mu_Ceres):

    from math import pi, sqrt
    mu_Sun = 132712E15 #m3/s2, nasa fact sheet
    a_Mer = 57.91E9 #m, nasa fact sheet
    e_Mer = 0.2056 #-, nasa fact sheet
    a_Ceres = 2.768*1.496E11 #m, nasa fact sheet * AU from google
    P_Mer = 87.97*24.0*3600.0 #s, google
    century = 100.0*365.25*24.0*3600.0 #s

    R_Mer_Ceres = a_Ceres - a_Mer
    orbitsPerCentury = century/P_Mer

    precessionPerOrbit = 1.5*pi \
        * ( mu_Ceres / mu_Sun ) \
        * ( (a_Mer / R_Mer_Ceres)**3.0 ) \
        * sqrt( 1.0 - e_Mer )

    precessionPerCentury = orbitsPerCentury * precessionPerOrbit

    return precessionPerCentury


from math import pi
a_Mer = 57.91E9 #m, nasa fact sheet
circ_merc = pi*a_Mer
radToArcsec = 3600.0*180.0/pi

p_J2_rad = precessionDueToJ2(2.25E-9)
p_J2_arcsec = p_J2_rad*radToArcsec

mu_Ceres = 62.78E9  # m3/s2, Kuchynka & Folkner 2013
p_Ceres_rad = precessionDueToCeres(mu_Ceres)
p_Ceres_arcsec = p_Ceres_rad*radToArcsec

print(p_J2_arcsec)
print(p_Ceres_arcsec)
print(precessionDueToCeres(mu_Ceres*1.01)*radToArcsec)




print("advance per year:", circ_merc*(p_J2_arcsec/3600.0)/100.0)

