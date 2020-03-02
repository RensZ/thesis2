"""
Created by Rens van der Zwaard on 2020-3-2

Purpose: a place to put self-made tools that are used in other files

"""


def Knm(degree,order):

    # reproduced function from TUDAT, as seen in line 761 of:
    # http://doxygen.tudat.tudelft.nl/d2/d5e/legendre_polynomials_8cpp_source.html

    from math import factorial, sqrt

    if order == 0:
        delta = 1.0
    else:
        delta = 0.0

    numerator = factorial(degree + order)
    denominator =  (2.0-delta) * (2*degree+1) * factorial(degree-order)
    norm = 1.0/sqrt(numerator/denominator)

    return norm


