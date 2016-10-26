"""
For the spectroscopic twins from Maedler et al (arXiv:1606.03015) calculate various quantities to test
whether the magnitude differences are consistent with the parallax differences.

Idea from the "Not ready for Gaia" session at the Gaia Sprint NYC.

THIS VERSION uses 2MASS photometry in addition to the Gaia G band, in order to account for exinction
(following Maedler et al).

Anthony Brown Oct 2016
"""

from __future__ import print_function

import numpy as np
import argparse
import os

from astropy.io.votable import parse
from astropy.table import Table, Column

def calculate_output_row(pleiad, table, nameA, nameB, dataA, dataB):
    """
    Calculate one row in the output table, corresponding to a pair of twins for a given Pleiades member.

    Parameters
    ----------

    pleiad - ID of Pleiades member.
    table - Table with results, new row will be appended after the calculations.
    nameA - ID of component A of pair of stars.
    nameB - ID of component B of pair of stars.
    dataA - Data for component A of pair (row from input VO table).
    dataB - Data for component B of pair (row from input VO table).
    """
    fluxAParallaxBSqr = dataA['phot_g_mean_flux']*dataB['parallax']**2
    fluxBParallaxASqr = dataB['phot_g_mean_flux']*dataA['parallax']**2

    errAB = fluxAParallaxBSqr * np.sqrt( (dataA['phot_g_mean_flux_error']/dataA['phot_g_mean_flux'])**2 +
            4*(dataB['parallax_error']/dataB['parallax'])**2)
    errBA = fluxBParallaxASqr * np.sqrt( (dataB['phot_g_mean_flux_error']/dataB['phot_g_mean_flux'])**2 +
            4*(dataA['parallax_error']/dataA['parallax'])**2)
    
    delta_lambda = (fluxAParallaxBSqr-fluxBParallaxASqr)/np.sqrt(errAB**2+errBA**2)
    
    deltaG = dataA['phot_g_mean_mag'] - dataB['phot_g_mean_mag']
    deltaDM = 5*np.log10(dataB['parallax']) - 5*np.log10(dataA['parallax'])

    gErrA = 2.5/np.log(10)*dataA['phot_g_mean_flux_error']/dataA['phot_g_mean_flux']
    gErrB = 2.5/np.log(10)*dataB['phot_g_mean_flux_error']/dataB['phot_g_mean_flux']
    deltaG_err = np.sqrt(gErrA*gErrA+gErrB*gErrB)
    
    errLogPlxA = 5/np.log(10)*dataA['parallax_error']/dataA['parallax']
    errLogPlxB = 5/np.log(10)*dataB['parallax_error']/dataB['parallax']
    deltaDM_err = np.sqrt(errLogPlxA*errLogPlxA+errLogPlxB*errLogPlxB)
    
    deltaH = dataA['h_m'] - dataB['h_m'] - 0.3 * ( (dataA['h_m']-dataA['ks_m']) - \
            (dataB['h_m']-dataB['ks_m']) )
    change_flux_A = dataA['phot_g_mean_flux'] - dataB['phot_g_mean_flux'] * \
            (dataA['parallax']**2/dataB['parallax']**2)

    change_parallax_A = dataA['parallax'] \
            - dataB['parallax'] * np.sqrt(dataA['phot_g_mean_flux']/dataB['phot_g_mean_flux'])

    change_parallax_A_tmass = dataA['parallax'] - dataB['parallax']*10**(-0.2*deltaH)

    table.add_row( (pleiad, nameA, nameB, fluxAParallaxBSqr, fluxBParallaxASqr, errAB, errBA,
        delta_lambda, deltaG, deltaDM, deltaG_err, deltaDM_err, deltaH, change_flux_A,
        dataA['phot_g_mean_flux_error'], change_parallax_A, change_parallax_A_tmass, dataA['parallax'],
        dataA['parallax_error'], dataB['parallax'], dataB['parallax_error']) )

def process_data(args):
    """
    For each spectroscopic pair (A, B) calculate flux_A*parallax_B^2 vs flux_B*parallax_A^2. These values
    should be (are assumed to be) the same for spectroscopic twins. In adddition calculate the magnitude
    difference and its error, as well as the difference in 5*log10(parallax) and its error.

    Pleiades stars, their candidate twins, and their identifications are taken from tables 1 and 2 in
    Maedler et al. For each pair among the pleiades member and its twin candidates the above quantities
    are calculated (thus including candidate-candidate pairing).

    Parameters
    ----------
    
    args - command line arguments

    Outputs
    -------

    Fits table with the following columns:

    Pleiades_ID : Name of Pleiades star for which candidate spectroscopic twins are listed in table 1 of
    Maedler et al
    ID_A : ID of component A of a pair of twins
    ID_B : ID of component B of a pair of twins
    lambdaAB : Value of flux_A*parallax_B^2
    lambdaBA : Value of flux_B*parallax_A^2
    lambdaAB_err : Error on flux_A*parallax_B^2
    lambdaBA_err : Error on flux_B*parallax_A^2
    delta_lambda : Error normalized difference of lambdaAB and lambdaBA
    deltaG : Value of G_A-G_B (EXTINCTION NOT ACCOUNTED FOR)
    deltaDM : Value of 5*log10(parallax_B)-5*log10(parallax_A)
    deltaG_err : Error on deltaG
    deltaDM_err : Error on deltaDM
    deltaH : Value of EXTINCTION CORRECTED difference in H-band magnitude
    change_flux_A : By how much does the flux of component A need to be changed if we assume the parallax
    ratio is correct?
    flux_A_err : The error on the flux of component A
    change_parallax_A : By how much does the parallax of component A need to be changed if we assume the
    flux ratio is correct?
    change_parallax_A_tmass : By how much does the parallax of component A need to be changed if we assume the
    2MASS extinction corrected magnitude difference to be correct?
    parallax_A : Parallax of component A
    parallax_A_err : Error on the parallax of component A
    parallax_B : Parallax of component B
    parallax_B_err : Error on the parallax of component B
    """

    # Pleiades star with Hipparcos identifiers.
    twinship = parse("data/pleiades-twins-hip-2mass.vot",
            pedantic=False).get_first_table().to_table(use_names_over_ids=True)
    # Pleiades stars with HII identifiers.
    twinshii = parse("data/pleiades-twins-hii-2mass.vot",
            pedantic=False).get_first_table().to_table(use_names_over_ids=True)

    # Twin candidates for each of the Pleiades stars with a Hipparcos ID.
    twinshipdict = {'HIP16979' : [43299, 1481, 3924],
            'HIP17091': [43299, 36312, 2751, 3924, 2724],
            'HIP17044': [6572, 95149],
            'HIP17316': [1481, 1825],
            'HIP17289': [5099, 41282, 17838, 18658, 29932, 19877],
            'HIP17511': [22844, 1427, 37844, 46934],
            'HIP18955': [41282, 5709, 7443, 20350]}

    # Twin candidates for each of the Pleiades stars with an HII ID.
    twinshiidict = {'HII430' : [14684], 
            'HII1215' : [3203, 41587, 109110, 25002],
            'HII1593' : [116819, 48141, 38041, 91700, 107805],
            'HII1794' : [95149],
            'HII1924' : [490],
            'HII2311' : [14684],
            'HII2506' : [490, 1825],
            'HII3179' : [33212, 112117, 72134, 45685, 38765, 413, 5280, 26722, 4747, 5806, 23128, 108859,
                115803, 116106, 3466, 53094]}

    output_table = Table(names=('Pleiades_ID', 'ID_A', 'ID_B', 'lambdaAB', 'lambdaBA', 'lambdaAB_err',
        'lambdaBA_err', 'delta_lambda', 'deltaG', 'deltaDM', 'deltaG_err', 'deltaDM_err', 'deltaH',
        'change_flux_A', 'flux_A_err', 'change_parallax_A', 'change_parallax_A_tmass', 'parallax_A',
        'parallax_A_err', 'parallax_B', 'parallax_B_err'),
        dtype=('S9', 'S9', 'S9', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f',
            'f', 'f', 'f', 'f'), meta={'name': 'Pleiades and their spectroscopic twins'})

    for pleiad in sorted(twinshipdict.keys()):
        index = np.argwhere(twinship['hip']==int(pleiad.split('HIP')[1]))[0][0]
        idlist = []
        if twinship[index]['parallax']:
            idlist.append(int(pleiad.split('HIP')[1]))
        else:
            continue
        for elem in twinshipdict[pleiad]:
            idlist.append(elem)
        for i in range(len(idlist)-1):
            for j in range(i+1,len(idlist)):
                idA = idlist[i]
                idB = idlist[j]
                searchA = np.argwhere(twinship['hip']==idA)
                if (searchA.size != 0):
                    dataA = twinship[searchA[0][0]]
                else:
                    continue
                searchB = np.argwhere(twinship['hip']==idB)
                if (searchB.size != 0):
                    dataB = twinship[searchB[0][0]]
                else:
                    continue
                calculate_output_row(pleiad, output_table, idA, idB, dataA, dataB)

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
        for i in range(len(idlist)-1):
            for j in range(i+1,len(idlist)):
                idA = idlist[i]
                idB = idlist[j]
                if (type(idA)==str):
                    searchA = np.argwhere(twinshii['hii']==bytes(idA, 'ascii'))
                    tableA = twinshii
                else:
                    searchA = np.argwhere(twinship['hip']==idA)
                    tableA = twinship
                if (searchA.size != 0):
                    dataA = tableA[searchA[0][0]]
                else:
                    continue
                if (type(idB)==str):
                    searchB = np.argwhere(twinshii['hii']==bytes(idB, 'ascii'))
                    tableB = twinshii
                else:
                    searchB = np.argwhere(twinship['hip']==idB)
                    tableB = twinship
                if (searchB.size != 0):
                    dataB = tableB[searchB[0][0]]
                else:
                    continue
                if type(idA)==str:
                    nameA = idA
                else:
                    nameA = 'HIP'+str(idA)
                if type(idB)==str:
                    nameB = idB
                else:
                    nameB = 'HIP'+str(idB)
                calculate_output_row(pleiad, output_table, nameA, nameB, dataA, dataB)

    output_table.write('twintables/pleiadestwins-2mass.fits', format='fits')

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser(description="""Data processing for spectroscopic twins test
            suggested at GaiaSprint NYC. This is for Pleaides stars and their twins.""")
    args=vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    process_data(args)
