
import numpy as np
from scipy import stats

class CountsSampler(object):
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

    def __init__(self, beta):
        self.beta    = beta

    def sample(self, counts):
        self.counts = counts
        self.weights = np.zeros( counts.shape[0] )
        self._update_weights()
        return self.weights

    def _update_weights( self ):
        counts_per_state = np.array(self.counts.sum(axis=1)).flatten() + 10.**-8
        w = np.power( counts_per_state, self.beta-1.0 ) # weigths
        w /= w.sum()
        self.weights = np.array(w).flatten()

    def get_sampling_states( self, num_samples ):
        self.update_weights()
        N = self.counts.shape[0]                               # num states
        sampler = stats.rv_discrete(name='sampler', values=[ np.arange(N), self.weights ])
        return sampler.rvs( size=num_samples )
