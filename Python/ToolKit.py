"""
Created by Rens van der Zwaard on 2020-3-2

Purpose: a place to put self-made tools that are used in other files

"""


def Knm(degree,order):

    """
    a function to calculate the normalization factor, needed to retreive the unnormalized J2 values
    """

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


def format_spines(ax, i, no_arcs):

    """
    a function to format axes with scheurlijnen
    """

    if no_arcs == 1:
        return

    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    d = .015

    #leftmost plot
    if i == 1:
        ax.spines['right'].set_visible(False)

        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    #middle plots
    if (i != 1) and (i != no_arcs):
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.set_yticklabels([])

        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
        kwargs.update(transform=ax.transAxes)
        ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax.plot((-d, +d), (-d, +d), **kwargs)

    #rightmost plots
    if i == no_arcs:
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(True)
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.tick_right()

        kwargs.update(transform=ax.transAxes)
        ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax.plot((-d, +d), (-d, +d), **kwargs)

    return