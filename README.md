# SpectroscopicTwins

Python code for spectroscopic twins plot produced for Gaia Sprint NYC

Take the spectroscopic twins for Pleiades stars listed in tables 1 and 2 of Maedler et al (2016,
https://arxiv.org/abs/1606.03015) and test whether their magnitude differences are consistent with
their Gaia DR1 parallax differences, on the assumption that for spectroscopic twins the absolute
magnitudes in the Gaia G-band are the same (ignoring extinction effects).

The test is done by calculating for each pair of twins (A, B) the value of flux_A * parallax_B^2 vs.
flux_B * parallax_A^2. These values should be the same under the above assumption. These quantities
have the advantadge of being linear in flux and that they have comparable error-bars. In addition
the more straightforward to understand differences in G and 5 * log10(parallax) are calculated, as
well as the difference in H (extinction corrected, as in equation 4 of Maedler et al).

For each Pleiades star the values are calculated for the star paired with each of its twin
candidates. In addition each of the twin candidates are also paired and the quantities above
calculated.

Dependencies
------------

Numpy, Astropy, Matplotlib

Data files
----------

pleiades-twins.csv: CSV file with Hipparcos identifiers of all the PelsXX Pleiades members and the twin
candidates.

HII-identifiers.txt: Pleiades member HII identifier and Gaia source_id

pleiades-twins-hii-2mass.vot: Table extracted form Gaia DR1 Archive by cross-matching through
pleiades-twins.csv, includes 2MASS photometry

pleiades-twins-hip-2mass.vot: Table extracted form Gaia DR1 Archive by cross-matching
HII-identifiers.txt, includes 2MASS photometry

Code
----

process-gaiadr1-2mass-data.py: Process the input data above and store result in a FITS file (in
twintables/ folder).

fancyaxes.py, distinct_colours.py: Plotting utilities.
yvsx-colourcoded.py: Quick plot of Y vs X from the output file in twintables/.
