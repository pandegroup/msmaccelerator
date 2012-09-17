
"""
A small library to test adaptive sampling routines.

To Do:
-- Optimize for speed
"""

import sys
import multiprocessing

import numpy as np
import scipy
from scipy import stats, sparse
import matplotlib.pyplot as plt

from msmbuilder import MSMLib


# --- observable_functions -----------------------------------------------------


def total_states(T):
    return T.shape[0]


def second_eigenvalue(T):
    v = np.real( sparse.linalg.eigs(T, k=5, which='LR', return_eigenvectors=False) )
    return v[1]


# --- functions for simulating adaptive sampling rounds ------------------------


def test_sampling_method(SamplerObject, observable_function=None,
                         transition_matrix=None, num_trials=5, rounds_of_sampling=50,
                         num_states=100, simultaneous_samplers=10,
                         size_of_intial_data=10, length_of_sampling_trajs=10,
                         procs=1):
    """
    Performs a test on the convergence of an adaptive sampling method, based on
    sampling from either a given or randomly generated "test" transition matrix.

    This mini-algorithm takes in a `SamplerObject`, per the MSMAccelerator
    conventions, and uses it to generate a simulation of many rounds of adaptive
    sampling. Sampled is a Markov chain, either user-specified (`transition_matrix`)
    or randomly generated. Many options for exactly how to perform sampling are
    provided.

    Finally, support is included for an `observable_function`, a function who's
    convergence properties will be calculated. This function should map a
    transition matrix to a scalar, and will typically be some value that the
    adaptive sampling routine is attempting to optimize.

    Parameters
    ----------
    SamplerObject : instantiated class
        An instance of a Sampler Object from MSMAccelerator's Sampling library

    observable_function : function
        A function mapping a transition matrix to a scalar, typically some value
        who's convergence we care about.

    transition_matrix : matrix
        The "true" matrix that will be used to generate simulated sampling. Pass
        None to generate a random matrix.

    procs : int
        Number of processors to utilize in parallel. Parallelism is over `trials`.


    Returns
    -------
    avg_distance : array_like, float
        The norm between the sampled transition matrix and the true
        (`transition_matrix`) matrix, averaged over all `num_trials`

    std_distance : array_like, float
        The standard deviation of the above norm.

    obs_distance_avg : array_like, float
        The difference between the value of `observable_function` for the
        sampled transition matrix and the true (`transition_matrix`) matrix,
        averaged over all `num_trials`

    obs_distance_std : array_like, float
        The standard deviation of the above difference.
    

    Other Parameters
    ----------------
    num_trials : int
        Tells the algorithm to perform this number of repeat sampling simulations.
        Each "trial" starts fresh with zero observed counds, and performs the requested
        `rounds_of_sampling`.

    rounds_of_sampling : int
        The number of rounds of adaptive sampling to do for each "trial".

    num_states : int
         If 'transition_matrix' == None, then this indicates how big the random matix
         should be.

    simultaneous_samplers : int
         The number of independent samplers to run for each `rounds_of_sampling`. 

    size_of_initial_data : int
         How much sampling should be collected before beginning adaptive rounds.

    length_of_sampling_trajs : int
         How long adaptively sampled trajectories should be (in units of lag time).
    """

    # print out some useful info
    info_str = \
"""
 --- Adaptive Sampling Simulation ---
Trials :                         %d
Rounds of sampling per trial :   %d
Processors utilized :            %d

Calculating... (this may take a while)
""" % (num_trials, rounds_of_sampling, procs)
    print info_str

    # package relevant args into dict - args for each trial are identical
    arg_dict = {
        'SamplerObject' : SamplerObject,
        'observable_function' : observable_function,
        'transition_matrix' : transition_matrix,
        'rounds_of_sampling' : rounds_of_sampling,
        'num_states' : num_states,
        'simultaneous_samplers' : simultaneous_samplers,
        'size_of_intial_data' : size_of_intial_data,
        'length_of_sampling_trajs' : length_of_sampling_trajs,
        }

    # choose whether or not to use multiprocessing based on procs arg.
    if procs > 1:
        pool = multiprocessing.Pool( processes=min(procs, num_trials) )
        multi_args = (arg_dict,) * num_trials
        result = pool.map_async(_run_trial, multi_args)
        result.wait()
        ret = np.array(result.get())

        if observable_function:
            avg_distance = ret.mean(0)[0,:]
            std_distance = ret.std(0)[0,:]
            obs_distance_avg = ret.mean(0)[1,:]
            obs_distance_std = ret.std(0)[1,:]
        else:
            avg_distance = ret.mean(0)
            std_distance = ret.std(0)

    else: # no parallelism
        distance_to_target = np.zeros((num_trials, rounds_of_sampling))
        if observable_function:
            obs_distance = np.zeros((num_trials, rounds_of_sampling))
            
        for trial in range(num_trials):
            print " --- Trial %d --- " % (trial+1,)
            distance_to_target[trial,:], obs_distance[trial,:] = _run_trial(arg_dict)
        
        # after all is done, average the results over trials
        avg_distance = distance_to_target.mean(axis=0)
        std_distance = distance_to_target.std(axis=0)

        if observable_function:
            obs_distance_avg = obs_distance.mean(axis=0)
            obs_distance_std = obs_distance.std(axis=0)

    # report our results!
    if observable_function:
        return avg_distance, std_distance, obs_distance_avg, obs_distance_std
    else:
        return avg_distance, std_distance


def _run_trial(arg_dict):

    # inject the arg_dict into the local namespace - may be a bad idea...
    for key in arg_dict.keys():
        exec(key + " = arg_dict['" + key + "']")

    # initialize data structures to hold output
    distance_to_target = np.zeros(rounds_of_sampling)
    obs_distance = np.zeros(rounds_of_sampling)

    # the assignments array will hold all of the output of all simulations
    assignments = -1.0 * np.ones((rounds_of_sampling * simultaneous_samplers + 1,
                                      max(size_of_intial_data, length_of_sampling_trajs+1) ))

    # initialize the "true" transition matrix
    if not transition_matrix:
        assert num_states > 0
        C_rand = np.random.randint( 0, 100, (num_states, num_states) )
        C_rand += C_rand.T
        T = MSMLib.estimate_transition_matrix( C_rand )
    else:
        T = transition_matrix
        num_states = T.shape[0]
    T = sparse.csr_matrix(T)
    MSMLib.check_transition(T)
        
    if observable_function:
        try:
            obs_goal = observable_function(T)
        except Exception as e:
            print >> sys.stderr, e
            raise Exception("Error evaluating function: %s" % observable_function.__name__)
            
    assignments[0,:size_of_intial_data] = MSMLib.sample(T, None, size_of_intial_data)

    # iterate, adding simulation time
    for sampling_round in range(rounds_of_sampling):
        
        # apply the adaptive sampling method - we need to be true to what a
        # real simulation would actually see for the counts matrix
        mod_assignments = assignments.copy()
        mapping = MSMLib.renumber_states( mod_assignments )
        C_mod = MSMLib.get_count_matrix_from_assignments( mod_assignments )
        T_mod = MSMLib.estimate_transition_matrix(C_mod)
        adaptive_sampling_multivariate = SamplerObject.sample(C_mod)

        # choose the states to sample from (in the original indexing)
        state_inds = np.arange(len(adaptive_sampling_multivariate))
        sampler = stats.rv_discrete(name='sampler', 
                                    values=[state_inds, adaptive_sampling_multivariate])
        starting_states = sampler.rvs( size=simultaneous_samplers )
        starting_states = mapping[starting_states]

        # start new 'simulations' in each of those states
        for i,init_state in enumerate(starting_states):
            a_ind = sampling_round * simultaneous_samplers + i + 1
            s_ind = length_of_sampling_trajs + 1
            assignments[a_ind,:s_ind] = MSMLib.sample(T, init_state, s_ind)

        # build a new MSM from all the simulation so far
        C_raw = MSMLib.get_count_matrix_from_assignments( assignments, NumStates=num_states )
        C_raw = C_raw + C_raw.T # might want to add trimming, etc.
        T_pred = MSMLib.estimate_transition_matrix(C_raw) 

        # calculate the error between the real transition matrix and our best prediction
        assert T.shape == T_pred.shape
        distance_to_target[sampling_round] = np.sqrt( ((T_pred - T).data ** 2).sum() ) \
                                             / float(num_states)

        if observable_function:
            obs_distance[sampling_round] = np.abs(observable_function(T_mod) - obs_goal)

    return distance_to_target, obs_distance


def plot_sampling_results(avg_distance, std_distance,
                          obs_distance_avg=None, obs_distance_std=None, 
                          obs_name=None):


    fig = plt.figure()
    ax1 = plt.subplot(111)
    x = np.arange( len(avg_distance), dtype=int )

    ax1.errorbar( x, avg_distance, yerr=std_distance, color='b' )
    ax1.plot( x, avg_distance, lw=2, color='b' )

    ax1.set_xlabel('Adaptive Sampling Round')
    ax1.set_ylabel('Error in Transition Matrix')

    if obs_distance_avg != None:
        ax2 = plt.twinx()
        ax2.errorbar( x, obs_distance_avg, yerr=obs_distance_std, color='g' )
        ax2.plot( x, obs_distance_avg, lw=2, color='g' )
        ax2.set_ylabel('Error in %s' % obs_name)

    ax1.set_yscale('log')
    ax2.set_yscale('log')

    if obs_name:
        fn = '%s_convergence.pdf' % obs_name
    else:
        fn = 'sampling_convergence.pdf'

    plt.savefig(fn)
    print "Saved: %s" % fn
    plt.show()
    return


    
