"""
Plot Y vs X from data in twintables/pleiadestwins-2mass.fits. Points are colour coded by the Pleiades ID
in table 2 from Maedler et al.

Anthony Brown Oct 2016
"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import os
from cycler import cycler

from astropy.table import Table, Column
from fancyaxes import get_basic_xy_axis

color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896',
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
rc('axes', facecolor='ffffff')
rc('axes', prop_cycle=(cycler('color', color_sequence)))
rc('axes', axisbelow=True)
rc('xtick', direction='out')
rc('ytick', direction='out')
rc('grid', color='cbcbcb')
rc('grid', linestyle='-')
rc('grid', linewidth=0.5)
rc('grid', alpha=1.0)
rc('figure', facecolor='ffffff')
rc('figure', dpi=80)
rc('figure.subplot', bottom=0.125)

def make_plot(args):
    """
    For each spectroscopic pair (A, B) plot flux_A*parallax_B^2 vs flux_B*parallax_A^2. These values
    should be (are assumed to be) the same for spectroscopic twins, to within extinction effects.

    Parameters
    ----------
    
    args - command line arguments
    """
    twinsdata = Table.read('twintables/pleiadestwins-2mass.fits')

    # This dictionary to be completed (and moved to separate module)
    axislabels = {'lambdaAB' : '$\\lambda_\\mathrm{AB}=f_\\mathrm{A}\\varpi_\mathrm{B}^2$',
            'lambdaBA' : '$\\lambda_\\mathrm{BA}=f_\\mathrm{B}\\varpi_\mathrm{A}^2$',
            'deltaDM' : '$\\Delta(5\\log_{10}\\varpi)$',
            'deltaH' : '$\\Delta H$ (corrected for extinction)',
            'deltaG' : '$\\Delta G$ (NOT corrected for extinction)',}

    pleiads = np.unique(twinsdata['Pleiades_ID'])
    numpairs = []
    for pleiad in pleiads:
        numpairs.append(np.count_nonzero(twinsdata['Pleiades_ID']==pleiad))
    numpairs=np.array(numpairs)

    fig = plt.figure(figsize=(12,12))
    ax = get_basic_xy_axis()

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params('both', width=1.5, which='major', labelsize=18)

    if args['xlog']:
        ax.set_xscale("log")
    if args['ylog']:
        ax.set_yscale("log")

    x = twinsdata[args['abscissa']]
    y = twinsdata[args['ordinate']]
    x_err = np.zeros_like(x)
    y_err = np.zeros_like(y)
    if args['abscissaErr']:
        x_err = twinsdata[args['abscissaErr']]
    if args['ordinateErr']:
        y_err = twinsdata[args['ordinateErr']]

    for pleiad in pleiads[np.argsort(numpairs)[::-1]]:
        indices = np.where(twinsdata['Pleiades_ID']==pleiad)
        ax.errorbar(x[indices], y[indices], xerr=x_err[indices], yerr=y_err[indices], fmt='o',
                label=pleiad, mec='none')

    ax.legend(loc='best')

    ax.plot(ax.get_xlim(), ax.get_ylim(), 'k-', zorder=0)
    ax.set_xlabel(axislabels[args['abscissa']], fontsize=20)
    ax.set_ylabel(axislabels[args['ordinate']], fontsize=20)
    ax.set_title('Spectroscopic twins from M\\"adler et al. (arXiv:1606.03015)', fontsize=20)

    basename='twins'
    if args['pdfOutput']:
        plt.savefig(basename+'.pdf')
    elif args['pngOutput']:
        plt.savefig(basename+'.png')
    else:
        plt.show()

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser(description="""Summary plot for Pleiades spectroscopic twins.""")
    parser.add_argument("abscissa", type=str, help="""Table field along X axis""")
    parser.add_argument("ordinate", type=str, help="""Table field along Y axis""")
    parser.add_argument("--abscissaErr", dest="abscissaErr", required=False, type=str, help="""Error bar field along X axis""")
    parser.add_argument("--ordinateErr", dest="ordinateErr", required=False, type=str, help="""Error bar field along Y axis""")
    parser.add_argument("--xlog", action="store_true", dest="xlog", help="Logarithmic X-axis")
    parser.add_argument("--ylog", action="store_true", dest="ylog", help="Logarithmic Y-axis")
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
    args=vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    make_plot(args)
