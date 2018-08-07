neufit
------

neufit.py fits the expected long-term neutral community composition derived by Sloan et al. (2006) to species abundance tables, such as OTU abundance tables obtained by 16S sequencing of microbial communities.

neusim.py simulates the underlying stochastic death-birth process and can be used to generate neutral example data.

Both the fitting procedure and the simulation are used in Sieber et al. (2018).

Usage
~~~~~~~~~~

To fit the neutral prediction to the provided simulation data and join with the mock taxonomy run:

::

$ python neufit.py sim_data.csv -t sim_taxonomy.csv

To run a neutral simulation with 50 communities and 100 species for 100000 time steps:

::

$ python neusim.py

This will generate sample data every 1000 time steps. Use -h to display the simulation parameters and how to change them.

References
~~~~~~~~~~

Sloan, W. T., Lunn, M., Woodcock, S., Head, I. M., Nee, S. and Curtis, T. P. (2006). Quantifying the roles of immigration and chance in shaping prokaryote community structure. Environmental Microbiology, 8:732-740. https://doi.org/10.1111/j.1462-2920.2005.00956.x

Sieber, M., Pita, L., Weiland-Bräuer, N., Dirksen, P., Wang, J., Mortzfeld, B., Franzenburg, S., Schmitz, R. A., Baines, J. F., Fraune, S., Hentschel, U., Schulenburg, H., Bosch, T. C. G. and Traulsen, A. (2018). The Neutral Metaorganism. bioRxiv. https://doi.org/10.1101/367243
