"""
Test for the spectroscopic twins from Maedler et al (arXiv:1606.03015) whether the magnitudes differences
are consistent with the parallax differences.

Idea from the "Not ready for Gaia" session at the Gaia Sprint NYC.

Anthony Brown Oct 2016
"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import os
from cycler import cycler

from astropy.io.votable import parse
from astropy.table import Table, Column

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
    should be (are assumed to be) the same for spectroscopic twins.

    Pleiades stars, their twins, and their identifications are taken from tables 1 and 2 in Maedler et
    al. For each pair among the pleiades member and its twin candidates the values of flux*parallax^2 are
    calculated (thus including candidate-candidate pairing).

    Parameters
    ----------
    
    args - command line arguments
    """
    twinship = parse("pleiades-twins-hip.vot",
            pedantic=False).get_first_table().to_table(use_names_over_ids=True)
    twinshii = parse("pleiades-twins-hii.vot",
            pedantic=False).get_first_table().to_table(use_names_over_ids=True)

    twinshipdict = {16979 : [43299, 1481, 3924],
            17091: [43299, 36312, 2751, 3924, 2724],
            17044: [6572, 95149],
            17316: [1481, 1825],
            17289: [5099, 41282, 17838, 18658, 29932, 19877],
            17511: [22844, 1427, 37844, 46934],
            18955: [41282, 5709, 7443, 20350]}

    twinshiidict = {'HII430' : [14684], 
            'HII1215' : [3203, 41587, 109110, 25002],
            'HII1593' : [116819, 48141, 38041, 91700, 107805],
            'HII1794' : [95149],
            'HII1924' : [490],
            'HII2311' : [14684],
            'HII2506' : [490, 1825],
            'HII3179' : [33212, 112117, 72134, 45685, 38765, 413, 5280, 26722, 4747, 5806, 23128, 108859,
                115803, 116106, 3466, 53094]}

    resultship = {}
    for pleiad in sorted(twinshipdict.keys()):
        index = np.argwhere(twinship['hip']==pleiad)[0][0]
        idlist = []
        if twinship[index]['parallax']:
            idlist.append(pleiad)
        for elem in twinshipdict[pleiad]:
            idlist.append(elem)
        resultship[pleiad] = Table(names=('Pleiades_ID', 'ID_A', 'ID_B', 'flxAPlxBSqr', 'flxBPlxASqr',
            'flxAPlxBSqrErr', 'flxBPlxASqrErr', 'Delta'), dtype=('S7', 'S7', 'S7', 'f', 'f', 'f', 'f',
                'f'), meta={'name': pleiad})
        for i in range(len(idlist)-1):
            for j in range(i+1,len(idlist)):
                idA = idlist[i]
                idB = idlist[j]
                searchA = np.argwhere(twinship['hip']==idA)
                if (searchA.size != 0):
                    indexA = searchA[0][0]
                else:
                    continue
                searchB = np.argwhere(twinship['hip']==idB)
                if (searchB.size != 0):
                    indexB = searchB[0][0]
                else:
                    continue
                fluxAParallaxBSqr = twinship['phot_g_mean_flux'][indexA]*twinship['parallax'][indexB]**2
                fluxBParallaxASqr = twinship['phot_g_mean_flux'][indexB]*twinship['parallax'][indexA]**2
                errAB = fluxAParallaxBSqr * np.sqrt(
                        (twinship['phot_g_mean_flux_error'][indexA]/twinship['phot_g_mean_flux'][indexA])**2
                        + 4*(twinship['parallax_error'][indexB]/twinship['parallax'][indexB])**2)
                errBA = fluxBParallaxASqr * np.sqrt(
                        (twinship['phot_g_mean_flux_error'][indexB]/twinship['phot_g_mean_flux'][indexB])**2
                        + 4*(twinship['parallax_error'][indexA]/twinship['parallax'][indexA])**2)
                resultship[pleiad].add_row( (pleiad, idA, idB, fluxAParallaxBSqr, fluxBParallaxASqr,
                    errAB, errBA, (fluxAParallaxBSqr-fluxBParallaxASqr)/np.sqrt(errAB**2+errBA**2)) )

    resultshii = {}
    for pleiad in sorted(twinshiidict.keys()):
        indexlist = np.argwhere(twinshii['hii']==bytes(pleiad, 'ascii'))
        if (indexlist.size > 0):
            index = indexlist[0][0]
        else:
            continue
        idlist = []
        if twinshii[index]['parallax']:
            idlist.append(pleiad)
        else:
            continue
        for elem in twinshiidict[pleiad]:
            idlist.append(elem)
        resultshii[pleiad] = Table(names=('Pleiades_ID', 'ID_A', 'ID_B', 'flxAPlxBSqr', 'flxBPlxASqr',
            'flxAPlxBSqrErr', 'flxBPlxASqrErr', 'Delta'), dtype=('S7', 'S7', 'S7', 'f', 'f', 'f', 'f',
                'f'), meta={'name': pleiad})
        for i in range(len(idlist)-1):
            for j in range(i+1,len(idlist)):
                idA = idlist[i]
                idB = idlist[j]
                if (type(idA)==str):
                    searchA = np.argwhere(twinshii['hii']==bytes(idA, 'ascii'))
                else:
                    searchA = np.argwhere(twinship['hip']==idA)
                if (searchA.size != 0):
                    indexA = searchA[0][0]
                else:
                    continue
                if (type(idB)==str):
                    searchB = np.argwhere(twinshii['hii']==bytes(idB, 'ascii'))
                else:
                    searchB = np.argwhere(twinship['hip']==idB)
                if (searchB.size != 0):
                    indexB = searchB[0][0]
                else:
                    continue
                fluxAParallaxBSqr = twinship['phot_g_mean_flux'][indexA]*twinship['parallax'][indexB]**2
                fluxBParallaxASqr = twinship['phot_g_mean_flux'][indexB]*twinship['parallax'][indexA]**2
                errAB = fluxAParallaxBSqr * np.sqrt(
                        (twinship['phot_g_mean_flux_error'][indexA]/twinship['phot_g_mean_flux'][indexA])**2
                        + 4*(twinship['parallax_error'][indexB]/twinship['parallax'][indexB])**2)
                errBA = fluxBParallaxASqr * np.sqrt(
                        (twinship['phot_g_mean_flux_error'][indexB]/twinship['phot_g_mean_flux'][indexB])**2
                        + 4*(twinship['parallax_error'][indexA]/twinship['parallax'][indexA])**2)
                resultshii[pleiad].add_row( (pleiad, idA, idB, fluxAParallaxBSqr, fluxBParallaxASqr,
                    errAB, errBA, (fluxAParallaxBSqr-fluxBParallaxASqr)/np.sqrt(errAB**2+errBA**2)) )

    fig = plt.figure(figsize=(12,12))
    ax = plt.gca()
    ax.grid(True)
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params('both', width=1.5, which='major', labelsize=18)

    ax.set_xscale("log")
    ax.set_yscale("log")

    for pleiad in sorted(resultship.keys()):
        fluxAplxB = resultship[pleiad]['flxAPlxBSqr']
        fluxBplxA = resultship[pleiad]['flxBPlxASqr']
        fluxAplxBErr = resultship[pleiad]['flxAPlxBSqrErr']
        fluxBplxAErr = resultship[pleiad]['flxBPlxASqrErr']
        ax.errorbar(fluxAplxB, fluxBplxA, xerr=fluxAplxBErr, yerr=fluxBplxAErr, fmt='o', label=pleiad,
                mec='none')
        resultship[pleiad].write(str(pleiad)+'.dat', format='ascii.fixed_width')

    for pleiad in sorted(resultshii.keys()):
        fluxAplxB = resultshii[pleiad]['flxAPlxBSqr']
        fluxBplxA = resultshii[pleiad]['flxBPlxASqr']
        fluxAplxBErr = resultshii[pleiad]['flxAPlxBSqrErr']
        fluxBplxAErr = resultshii[pleiad]['flxBPlxASqrErr']
        ax.errorbar(fluxAplxB, fluxBplxA, xerr=fluxAplxBErr, yerr=fluxBplxAErr, fmt='o',
                label=pleiad, mec='none')
        resultshii[pleiad].write(pleiad+'.dat', format='ascii.fixed_width')

    ax.legend(loc='best')

    ax.plot(ax.get_xlim(), ax.get_ylim(), 'k-', zorder=0)
    ax.set_xlabel('$f_\\mathrm{A}\\varpi_\mathrm{B}^2$', fontsize=20)
    ax.set_ylabel('$f_\\mathrm{B}\\varpi_\mathrm{A}^2$', fontsize=20)
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
    parser = argparse.ArgumentParser(description="""Spectroscopic twins test suggested at GaiaSprint NYC.""")
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
    args=vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    make_plot(args)
