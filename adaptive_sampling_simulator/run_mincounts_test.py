from scipy import io
from msmaccelerator import Sampling
from TestAdaptive import *

T = io.mmread('ww_tProb.mtx')
beta = 0.0 # mincounts

counts_sampler = Sampling.counts_sampler(beta)

output = test_sampling_method(counts_sampler, observable_function=total_states, transition_matrix=T, procs=12)
plot_sampling_results( *output, obs_name='Number of States')
