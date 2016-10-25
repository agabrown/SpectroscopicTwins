# SpectroscopicTwins
Python code for spectroscopic twins plot produced for Gaia Sprint NYC

Take the spectroscopic twins for Pleiades stars listed in tables 1 and 2 of Maedler et al (2016,
https://arxiv.org/abs/1606.03015) and test whether their magnitude differences are consistent with their
Gaia DR1 parallax differences, on the assumption that for spectroscopic twins the absolute magnitudes in
the Gaia G-band are the same.

The test is done by plotting for each pair of twins (A, B) the value of flux_A * parallax_B^2 vs.
flux_B * parallax_A^2. These values should be the same under the above assumption. Plotting these ensures
that the error-bars in both axes have comparable values.

For each Pleiades star the values are calculated for the star paired with each of its twin candidates. In
addition each of the twin candidates are also paired and the product of flux and parallax^2 calculated.

Data files
----------

pleiades-twins.csv: CSV file with Hipparcos identifiers of all the PelsXX Pleiades members and the twin
candidates.

HII-identifiers.txt: Pleiades member HII identifier and Gaia source_id

pleiades-twins-hii.vot: Table extracted form Gaia DR1 Archive by cross-matching through
pleiades-twins.csv.

pleiades-twins-hip.vot: Table extracted form Gaia DR1 Archive by cross-matching HII-identifiers.txt.
