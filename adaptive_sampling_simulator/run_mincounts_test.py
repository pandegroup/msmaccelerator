from scipy import io
from msmaccelerator.adaptive import CountsSampler
from test_adaptive import *

T = io.mmread('ww_tProb.mtx')
beta = 0.0 # mincounts

counts_sampler = CountsSampler(beta) 

output = test_sampling_method(counts_sampler, observable_function=total_states, transition_matrix=T, procs=12)
plot_sampling_results( *output, obs_name='Number of States')
