# neusim: Simulate a neutral death-birth process
#
# The underlying stochastic process is described in 
# Sloan et al, Environ Microbiol 2006 8:732-740.
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
from argparse import ArgumentParser, ArgumentTypeError
from math import log10
from matplotlib import pyplot
from pandas import DataFrame
from random import random

def non_negative_int(arg):
    # Argparser type: non-negative int
    nnint = int(arg)
    if nnint < 0:
        raise ArgumentTypeError(arg + ' < 0, must be non-negative')
    return nnint
    
def float_0_1(arg):
    # Argparser type: float in the interval [0,1]
    f01 = float(arg)
    if f01 < 0 or f01 > 1:
        raise ArgumentTypeError(arg + ', m must be between 0 and 1')
    return f01

def random_species(proportions):
    # Pick a number (=species) based on weights (=species proportions in community)
    distance = random()
    for i, p in enumerate(proportions):
        distance -= p
        if distance < 0:
            return i

def colonize(source_community, n_community_size):
    # Build the initial communities

    n_species = source_community.size
    init_community = np.zeros(n_species, dtype=int)

    # Randomly pick n initial species (with replacement, so some might occur multiple times)
    init_species = np.random.choice(n_species, 1 + int(random()*(n_species-1)))

    # Divide community into equally sized chunks
    init_prop = n_community_size/len(init_species)
    for species in init_species:
        if n_community_size - init_community.sum() >= init_prop:
            init_community[species] += init_prop
    
    # Fill leftover empty space with last species, if any
    init_community[init_species[-1]] += n_community_size - init_community.sum()

    return init_community
    
# Parse command line arguments
parser = ArgumentParser(description='Simulate neutral death-birth process')
parser.add_argument('-t', dest='t_max', type=non_negative_int, default=100000, help='Simulation run time. Default: 100000')
parser.add_argument('-p', dest='t_sample', type=non_negative_int, default=1000, help='Sample interval. Default: 1000')
parser.add_argument('-m', dest='m', type=float_0_1, default=0.05, help='Immigration probability. Default: 0.05')
parser.add_argument('-c', dest='n_communities', type=non_negative_int, default=50, help='Number of local communities. Default: 50')
parser.add_argument('-s', dest='n_community_size', type=non_negative_int, default=1000, help='Local community size. Default: 1000')
parser.add_argument('-o', dest='n_species', type=non_negative_int, default=100, help='Number of species. Default: 100')
args = parser.parse_args()

t_max = args.t_max
t_sample = args.t_sample
m = args.m
n_communities = args.n_communities
n_community_size = args.n_community_size
n_species = args.n_species


# Define the source community (e.g. here drawn from a geometric distribution)
source_community = np.random.geometric(0.00001, n_species)
source_community = 1.0*source_community/source_community.sum()

# Initialize the local communities
local_communities = np.array([colonize(source_community, n_community_size) for i in range(n_communities)])

# Run the simulation for tmax time steps
for t in range(t_max+1):

    # One local death-birth event per community
    for community in local_communities:
        # Death
        community[random_species(1.0*community/n_community_size)] -= 1
        # Birth
        community[random_species(m*source_community + (1-m)*community/(n_community_size-1))] += 1

    # Save current communities
    if t % t_sample == 0:
        sample_id = ['sample_' + str(i+1) for i in range(len(local_communities))]
        species_id = ['otu_' + str(i+1) for i in range(len(source_community))]
        DataFrame(np.transpose(local_communities), columns=sample_id, index=species_id).to_csv('abundances_t' + str(t).zfill(int(log10(t_max))+1)  + '.csv', sep='\t')
        print 'sampled at t=' + str(t)
