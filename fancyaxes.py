"""
Provides styled axes (plot boxes) inspired on the R default plotting style.

Anthony Brown Aug 2015
"""

from __future__ import print_function

import matplotlib.pyplot as plt
from matplotlib import rc, cm
from cycler import cycler

from distinct_colours import get_distinct

#colour_sequence = get_distinct(12)
colour_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896',
        '#9467bd', '#c5b0d5', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7',
        '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
#color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896',
#        '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7',
#        '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

# Configure matplotlib
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amsmath}')
rc('font', family='serif', size=16)
rc('xtick.major', size='6')
rc('xtick.minor', size='4')
rc('ytick.major', size='6')
rc('ytick.minor', size='4')
rc('lines', linewidth=1.5)
rc('axes', linewidth=1)
rc('axes', facecolor='f9f9f9')
rc('axes', prop_cycle=(cycler('color',colour_sequence)))
rc('xtick', direction='out')
rc('ytick', direction='out')
rc('grid', color='cbcbcb')
rc('grid', linestyle='-')
rc('grid', linewidth=0.5)
rc('grid', alpha=1.0)
rc('figure', facecolor='ffffff')
rc('figure', dpi=80)
rc('figure.subplot', bottom=0.125)

def get_basic_xy_axis(withgrid=True):
    """
    Obtain an Axes object for basic XY plots.

    Parameters
    ----------

    None

    Keywords
    --------

    withgrid - When true a grid is displayed in the plot background

    Returns
    -------

    Styled Axes object which can be used for further plotting instructions.
    """

    ax = plt.gca()
    configure_basic_xy_axis(ax, withgrid=withgrid)
    return ax

def configure_basic_xy_axis(ax, withgrid=True):
    """
    Apply a basic configuration to the input axis object.

    Parameters
    ----------

    ax - The axis object to configure.

    Keywords
    --------

    withgrid - When true a grid is displayed in the plot background

    Returns
    -------

    Nothing
    """

    if withgrid:
        ax.grid(True)
    # Move left and bottom spines outward by 10 points
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
