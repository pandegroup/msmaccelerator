
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

#def even(ff_data, num_jobs):
#    """ For n different forcefields, sample num_jobs / n from the uniform distrubtion
#    over the microstates for each forcefield """
#
#    output = {}
#    jobs_generated = 0 # number generated so far
#    num_forcefields = len(ff_data)
#
#    for i, ff in enumerate(ff_data):
#        num_jobs_this_ff = (num_jobs - jobs_generated) / (num_forcefields - i)
#        num_microstates = ff['counts'].shape[0]
#        jobs_generated += num_jobs_this_ff
#
#        output[ff['name']] = np.random.randint(0, num_microstates, size=num_jobs_this_ff)
#
#    return output

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
        sampler = counts_sampler(beta, ff['counts'])
        microstate_weights[i] = sampler.weights

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
# The big benefit of this method is that it allows for simple ways to combine 
# methods when one may want to optimize many features at once. A simple linear
# combination of weights yields a new categorical distribution that contains
# information from all included methods, allowing the workflow to proceed in
# a smooth and simple fashion.
#

class counts_sampler(object):
    """
    A method for counts-based sampling. Involves a temperature
    factor 'beta' that controls the level of exploration vs. 
    refinement. Works as follows:

      beta = 0 (high temp) :   Full exploration, put simulations
                               where few counts have been obs.
      beta = 1 (room temp) :   Random sampling, no counts preferences
      beta > 1 (low  temp) :   Refinement, such that we focus on
                               areas with a high number of counts

    The explicit formula used is:

         Prob( choose state i ) ~ \sum_j C_{ij} ^{ beta - 1 }

    """

    def __init__(self, beta, counts):
        self.beta    = beta
        self.counts  = counts
        self.weights = np.zeros( counts.shape[0] )
        self.update_weights()

    def update_counts( self, new_counts ):
        self.counts = new_counts
        self.update_weights()

    def update_beta(self, beta):
        self.beta = beta
        self.update_weights()

    def update_weights( self ):
        w = np.power( self.counts.sum(axis=0).flatten(), self.beta-1.0 ) # weigths
        w /= w.sum()
        self.weights = np.array(w).flatten()

    def get_sampling_states( self, num_samples ):
        self.update_weights()
        N = self.counts.shape[0]                               # num states
        sampler = stats.rv_discrete(name='sampler', values=[ np.arange(N), self.weights ])
        return sampler.rvs( size=num_samples )
    

class dirichlet(object):
    """ A collection of useful functions for dealing with Dirichlet distributions.
        While these functions were generated with a mind for applications to MSMs,
        they should be generally applicable. """

    def mv_beta( array ):
        return np.product( special.gamma( array ) ) / special.gamma( np.sum( array ) )

    def max_likelihood_params( counts ):
        """ Returns the most likely parameters given a set of counts
        """
        C_sums = np.asarray( counts.sum(axis=1) ).flatten()
        D = sparse.dia_matrix( (1.0 / C_sums ,0), C.shape).tocsr()
        T = D.dot(C)
        return T

    def max_posterior_params( counts ):
        """ Returns the parameters at the mode of the posterior

        Note: formula for the mode (x_i)

              x_i = ( counts_i + alpha_i - 1 ) / ( \sum_i ( counts_i + alpha_i ) - N )

        where N is the size of the counts matrix

        Arguments:
          counts - the observed counts
        """
       
        # check this for speed later 
        N = counts.shape[0]
        D = np.asarray( counts.sum(axis=1) ).flatten() + self.alpha - N
        T = (counts + self.alpha - 1.0) * (1.0 / D)
        return T

    
class telescopic_estimator(dirichlet):
    """
    Combines data from many forcefields to provide a best estimate
    of the transition matrix parameters for some gold standard
    forcefield.


    Argument Notes    

    uniform_prior:   sets the prior to be uniform. If false, then takes all
                     Dirichlet alpha values to be 1/N


    Important Latent Variables

    alpha :          the Dirichlet prior parameter, which is taken to be a
                     constant across all FFs 

    """


    def __init__(self, uniform_prior=False):
        self.num_ff    = 0         
        self.uniform_prior = uniform_prior
        self.ff_counts = {}        # dict to hold all counts data
        self.gold      = ''        # gold standard FF's name
        self.Z         = {}        # the coupling parameters
        self.Z_hist    = {}        # the history of the z parameters

    def add_ff( self, name, counts, is_gold=False ):

        if name in ff_counts.keys():
            print "WARNING: %s already in ff_counts for telescopic_estimator" % name
        else:
            self.num_ff += 1

        self.ff_counts[ name ] = counts

        if is_gold:
            self.gold = name


    def update_counts( self, name, counts ):
        assert name in ff_counts.keys()
        self.ff_counts[ name ] = counts

    def _initialize_alpha(self):
        if self.uniform_prior:
            self.alpha = 1.0
        else:
            self.alpha = 1.0 / self.ff_counts[ self.gold ].shape[0]

    def estimate_gold_counts( self ):

        assert self.num_ff > 0 
 
        self._initialize_alpha()
        C_hat = self.ff_counts[ self.gold ] # start with the obs. gold counts

        for ff in self.ff_counts.keys():
            self._estimate_z( name )
            C_hat += self.ff_counts[name] * self.Z[name]

        C = self.max_posterior_params( C_hat )
        return C

    def _add_Zhist( self, z, name ):
        if name in self.Z_hist.keys():
            self.Z_hist[name].append(z)
        else:
            self.Z_hist[name] = [z]

    def _estimate_z( self, name, verbose=False ):
        """ Provides the MLE estimate for the coupling parameters z_i
            for the specified forcefield 'name' """

        # first, maximize the likelihood of FF 'X' with the given prior
        T_X = self.max_posterior( self.ff_counts[name] )
        N_X = np.sum( self.ff_counts[name], axis=1)       # the num. observations 
        C_X = T_X * N_X                                   # estimated counts

        # next, optimize the log likelihood
        N     = counts.shape[0]
        z_opt = np.zeros( N )

        guess = 1.0 / C_X.sum()
        for i in range(N):

            C_X_i = C_X[i,:]
            C_Y_i = self.ff_counts[self.gold][i,:]

            def objective( z ):
                """ the log-likelihood objective function - minize me """
                LL = np.sum( special.gammaln( C_Y_i + z * C_X_i + 1 ) ) \
                             - special.gammaln( z * C_X_i.sum() )
                return -LL

            def D_objective( z ):
                """ The (semi)-analyical derivative of the objective function """
                dLL = np.sum( C_X_i * special.psi( C_Y_i + z * C_X_i + 1 )) \
                              - C_X_i.sum() * special.psi( z * C_X_i.sum() )
                return -dLL

            z_opt[i], info = optimize.minimize( objective, guess, full_output=True,
                                                method='BFGS', jac=D_objective )
            if verbose: print info
            guess = z_opt[i] # set the previous answer as the initial guess for the next state
            
        self.Z[name] = z
        self._add_Zhist( z, name )


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
