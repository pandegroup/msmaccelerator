
"""
Each adaptive sampling algorithm is implemented as a function with the
following signature.

Inputs:
    ff_data: a list of dicts for each forcefield. Each dict
    represents a forcefield and should contain:
        'name' - the name of the forcefield
        'counts' - a sparse matrix of the counts for the MSM constructed in
                   this ff
        'cost' - a float represnting the relative cost of sampling in this
                 forcefield with respect to the other forcefields
        'gold_standard' - boolean to represnt whether this forcefield is
                          the "gold standard". Should only be true for one
                          forcefield    
    num_jobs: The number of jobs to generate
    
Returns: A dict where the keys are the name of each forcefield, and the 
values are a list of the microstates in that forcefield to sample from
"""

### Imports -------------------------------------------------------------------
import numpy as np

from scipy import stats
from scipy import sparse
from scipy import optimize


### Adaptive Sampling Routines ------------------------------------------------
# 
# Different adaptive sampling routines will be placed here. By routine, we mean
# a function that is essentially a map that takes 'counts' data, collected in
# different forcefields and produces a list of states to sample from.
#
# Specifically, they should be of the form:
# 
#              routine( ff_data, num_jobs ) --> { states }
#
# ff_data:  dict of the form ff_data[ 'ff_name' ] = counts_sparse_mtx
# num_jobs: integer indicating the number of states to generate jobs for,
#           *in all forcefields* (i.e. the total across FFs)
# states:   dict of form states[ 'ff_name' ] = list_of_states_to_sample
#
# These 'routines' may combine many adaptive 'methods', which are constructed
# below as objects. For instance, one might want to sample with a combination
# of counts and eigenvalue based sampling... they should write a routine
# specific to their application constructed from these individual general 
# classes.
#
# Final Note: When done writing a new routine, add your method to the simple
# test at the bottom of this file
#


def even(forcefields):
    num_forcefields = len(forcefields)
    
    # put all the weight on the first forcefield
    ff_weights = np.zeros(num_forcefields)
    ff_weights[0] = 1
    
    microstate_weights = [None] * num_forcefields
    for i, ff in enumerate(forcefields):
        num_microstates = ff['counts'].shape[0]
        microstate_weights[i] = np.ones(num_microstates) / num_microstates
    
    return ff_weights, microstate_weights


def explorative_counts(forcefields):
    """ This is a simple counts-based adaptive sampling algorithm that
        maximizes the size of explored space """

    beta = 0.0 # setting the temp to 0 means full explorer mode
    num_forcefields = len(forcefields)

    # put all the weight on the first forcefield
    ff_weights = np.zeros(num_forcefields)
    ff_weights[0] = 1.0

    microstate_weights = [None] * num_forcefields

    for i, ff in enumerate(forcefields):
        num_microstates = ff['counts'].shape[0]
        sampler = counts_sampler(beta)
        microstate_weights[i] = sampler.sample(ff['counts'])

    return ff_weights, microstate_weights


### Adaptive Sampling Methods -------------------------------------------------
#
# This section contains the general classes used to construct adaptive sampling
# routines. While there is a lot of flexibility in writing or adding to these
# classes, we have adopted the following convention that should be very general.
#
# Each class below is designed to produce some set of "weights" for each state,
# in general these weights can be though of as a discrete categorial distribution
# across the MSM space. Then, each state will be sampled with probability
# proportional to its weight. The weights are set by different methods to
# optimize some MSM parameter of interest, whether total states explored, 
# uncertainty in a specific eigenmode, etc.
#
# These routines have been standardized to be used in the following two-step
# fashion:
# 
# (1) Initialize:       sampler = SamplingMethod( args_controlling_algorithm )
# (2) Request weights:  multinomial = sampler.sample( transition_counts_matrix )
#
# The big benefit of this method is that it allows for simple ways to combine 
# methods when one may want to optimize many features at once. A simple linear
# combination of weights yields a new categorical distribution that contains
# information from all included methods, allowing the workflow to proceed in
# a smooth and simple fashion.
#


### Testing -------------------------------------------------------------------    
# a simple test to make sure function calls are returning properly, each 
# routine should be included here (probably)

if __name__ == '__main__':
    print even([{'name': 1,   'counts': sparse.eye(10,10)},
                {'name': 2,   'counts': sparse.eye(10,10)},
                {'name': 'a', 'counts': sparse.eye(10,10)}])
    print explorative_counts([{'name': 1,   'counts': sparse.eye(10,10)},
                              {'name': 2,   'counts': sparse.eye(10,10)},
                              {'name': 'a', 'counts': sparse.eye(10,10)}])
