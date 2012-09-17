from scipy import io
from msmaccelerator.adaptive import SinghalSampler
from test_adaptive import *

T = io.mmread('ww_tProb.mtx')

indices = [1] # use only the second eigenvalue (lambda-2 sampling)
weights = [1.0]

ss = SinghalSampler(indices, weights, prior='Jefferys')

output = test_sampling_method( ss, observable_function=second_eigenvalue, transition_matrix=T, procs=12 )
plot_sampling_results( *output, obs_name='Number of States' )



