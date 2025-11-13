# neufit: Fit a neutral community model to species abundances, e.g. from an OTU table
#
# For the theory behind this see Sloan et al, Environ Microbiol 2006 8:732-740.
# To run on the example simulation data: python neufit.py sim_data.csv
# To link with the mock taxonomy use the -t sim_taxonomy.csv option
#
# Copyright (C) 2018 Michael Sieber (sieber.ecoevo.de)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import pandas as pd

from argparse import ArgumentParser, ArgumentTypeError, FileType
from lmfit import Parameters, Model, fit_report
from math import log10
from matplotlib import pyplot
from os.path import splitext
from scipy.stats import beta
from statsmodels.stats.proportion import proportion_confint


def beta_cdf(p, N, m):
    # Expected long term distribution under the neutral model (truncated cumulative beta-distribution)
    return beta.cdf(1.0, N*m*p, N*m*(1.0-p)) - beta.cdf(1.0/N, N*m*p, N*m*(1.0-p))

def subsample(counts, depth):
    # Subsamples counts to uniform depth, dropping all samples without enough depth
    for sample in counts:
        if counts[sample].sum() >= depth:
            flattened = np.repeat(np.arange(counts[sample].size), counts[sample])
            subsample = np.random.choice(flattened, depth, replace=False)
            counts[sample] = np.bincount(subsample, minlength=counts[sample].size)
        else:
            print('dropping sample ' + sample + ' with ' + str(counts[sample].sum()) + ' reads < ' + str(depth))
            counts = counts.drop(sample, axis=1)
    return counts

def non_negative_int(arg):
    # Argparser type: non-negative int
    nnint = int(arg)
    if nnint < 0:
        raise ArgumentTypeError(arg + ' < 0, must be non-negative')
    return nnint


# Parse command line arguments
parser = ArgumentParser(description='Fit neutral expectation to OTU abundance data')
parser.add_argument('data_file', type=FileType('r'), help='The OTU abundance table')
parser.add_argument('-t', dest='taxonomy_file', type=FileType('r'), help='Link OTUs to taxonomic information if available')
parser.add_argument('-r', dest='rarefaction_level', type=non_negative_int, default=0, help='Set rarefaction level. Default: highest possible uniform read depth.')
parser.add_argument('-i', dest='ignore_level', type=non_negative_int, default=0, help='Ignore OTUs below this abundance threshold. Default: 0')
args = parser.parse_args()

# Read data
print('reading ' + args.data_file.name)
abundances = pd.read_table(args.data_file, header=0, index_col=0, sep='\t').astype(int)
abundances = abundances[abundances.sum(1) > args.ignore_level]
print('dataset contains ' + str(abundances.shape[1]) + ' samples (sample_id, reads):')
print(abundances.sum(0))

# Determine uniform read depth
if args.rarefaction_level == 0 or args.rarefaction_level > max(abundances.sum(0)):
    args.rarefaction_level = min(abundances.sum(0))
    print('rarefying to highest possible uniform read depth',)
else:
    print ('rarefying to custom rarefaction level',)
print('(' + str(args.rarefaction_level) + ' reads per sample)')

# Optionally subsample the abundance table, unless all samples already have the required uniform read depth
if not all(n_reads == args.rarefaction_level for n_reads in abundances.sum(0)):
    abundances = subsample(abundances, args.rarefaction_level)
    abundances = abundances[abundances.sum(1) > 0]

# Dataset shape
n_otus, n_samples = abundances.shape
n_reads = args.rarefaction_level

print('fitting neutral expectation to dataset with ' + str(n_samples) + ' samples and ' + str(n_otus) + ' otus')

# Calculate mean relative abundances and occurrence frequencies
mean_relative_abundance = (1.0*abundances.sum(1))/n_reads/n_samples
occurrence_frequency = (1.0*np.count_nonzero(abundances, axis=1))/n_samples

occurr_freqs = pd.DataFrame(mean_relative_abundance, columns=['mean_abundance'])
occurr_freqs.index.name = 'otu_id'
occurr_freqs['occurrence'] = occurrence_frequency
occurr_freqs = occurr_freqs.sort_values(by=['mean_abundance'])

# Join with taxonomic information (optional)
if args.taxonomy_file != None:
    taxonomy = pd.read_table(args.taxonomy_file, header=0, index_col=0, sep='\t')
    taxonomy.index.name = 'otu_id'
    occurr_freqs = occurr_freqs.join(taxonomy)

# Fit the neutral model
params = Parameters()
params.add('N', value=n_reads, vary=False)
params.add('m', value=0.5, min=0.0, max=1.0)
beta_model = Model(beta_cdf)
beta_fit = beta_model.fit(occurr_freqs['occurrence'], params, p=occurr_freqs['mean_abundance'])

# Report fit statistics
r_square = 1.0 - np.sum(np.square(occurr_freqs['occurrence'] - beta_fit.best_fit))/np.sum(np.square(occurr_freqs['occurrence'] - np.mean(occurr_freqs['occurrence'])))
print(fit_report(beta_fit))
print('R^2 = ' + '{:1.2f}'.format(r_square))

# Adding the neutral prediction to results
occurr_freqs['predicted_occurrence'] = beta_fit.best_fit
occurr_freqs['lower_conf_int'], occurr_freqs['upper_conf_int'] = proportion_confint(occurr_freqs['predicted_occurrence']*n_samples, n_samples, alpha=0.05, method='wilson')

# Save non-neutral otus (here simply determined by lying outside the confidence intervals)
above = occurr_freqs[occurr_freqs['occurrence'] > occurr_freqs['upper_conf_int']]
below = occurr_freqs[occurr_freqs['occurrence'] < occurr_freqs['lower_conf_int']]
pd.concat((above, below)).to_csv(splitext(args.data_file.name)[0] + '_nonneutral.csv')

# Prepare results plot
pyplot.xlabel('Mean relative abundance across samples', fontsize=18)
pyplot.xscale('log')
x_range = np.logspace(log10(min(occurr_freqs['mean_abundance'])/10), 0, 1000)
pyplot.xlim(min(x_range), max(x_range))
pyplot.xticks(fontsize=16)
pyplot.ylabel('Occurrence frequency in samples', fontsize=18)
pyplot.ylim(-0.05, 1.05)
pyplot.yticks(fontsize=16)

# Plot data points
pyplot.plot(occurr_freqs['mean_abundance'], occurr_freqs['occurrence'], 'o', markersize=10, fillstyle='full', color='black')

# Plot best fit
pyplot.plot(x_range, beta_cdf(x_range, n_reads, beta_fit.best_values['m']), '-', lw=5, color='darkred')
lower, upper = proportion_confint(beta_cdf(x_range, n_reads, beta_fit.best_values['m'])*n_samples, n_samples, alpha=0.05, method='wilson')
pyplot.plot(x_range, lower, '--', lw=2, color='darkred')
pyplot.plot(x_range, upper, '--', lw=2, color='darkred')
pyplot.fill_between(x_range, lower, upper, color='lightgrey')

pyplot.text(0.05, 0.9, '$R^2 = ' + '{:1.2f}'.format(r_square) + '$', fontsize=16, transform=pyplot.gca().transAxes)

pyplot.tight_layout()
pyplot.savefig(splitext(args.data_file.name)[0] + '_neutral_fit.png')
