# This file is part of MSMAccelerator.
#
# Copyright 2011 Stanford University
#
# MSMAccelerator is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

"""
sampling.py: Adaptive sampling routines

Adaptive sampling has been split into two components, for enhanced
flexibility and modularity:

Methods  - Adaptive "methods" are mathematical models that predict which
           states it would be most beneficial to sample based on some
           objective function. Examples include:
           -- Maximizing the number of states found
           -- Minimizing the error in the transition matrix
           -- Minimizing the error in the stationary probabilities

Routines - Adaptive "routines" are composed of methods, and reflect
           complete algorithms designed to execute a certain kind of
           simulation. An example routine might be:

           (1) Maximize states found using a method designed to do so. 
               Continue until some convergence criteria is met.
           (2) Employ a method designed to minimize the error in the 
               transition matrix.
           
           A simpler, more concrete version of this might be:

           (1) Sample from the least visited state for 100 iterations
           (2) Sample from each state with equal probability for 100
               iterations

Practically, "methods" are implemented as classes. They are called like:

(1) Initialize:       sampler = SamplingMethod( args_controlling_algorithm )
(2) Request weights:  multinomial = sampler.sample( transition_counts_matrix )

where `multinomial` is a vector giving the probability you should sample from
each state.

"Routines" are functions. The format of a "routine" should be as follows:



This file houses the routines. Methods are in the "Adaptive" module.
"""


from __future__ import division

import numpy as np
from scipy import stats
from scipy import sparse
from scipy import optimize

import logging
import IPython as ip
from msmaccelerator.models import Trajectory, Forcefield, MarkovModel, MSMGroup
from sqlalchemy.sql import and_, or_, not_
from sqlalchemy import func

import adaptive

logger = logging.getLogger('MSMAccelerator.sampling')


def myfavorite(Session, msmgroup):
    """ An example """
    logger.info('Running myfavorite')

    def activation_response(x, k):
        """ Curve from [0,1] -> [0,1] """
        if k < 0:
            k = k / (1 - k)
        return x / (1 - k * (x - 1))
    
    # =============#
    k_factor = -100
    # =============#
    prev = Session.query(MSMGroup).get(msmgroup.id - 1)
    if prev is None:
        return default(Session, msmgroup)
    
    # number of new states discovered
    n_new = msmgroup.n_states - prev.n_states
    if n_new < 0:
        return default(Session, msmgroup)
    
    q = Session.query(Trajectory)
    new_trajectories = q.filter(Trajectory.msm_groups.contains(msmgroup)).filter(not_(
            Trajectory.msm_groups.contains(prev)))

    # sum of the number of steps in the new trajectories
    sum_op = Session.query(func.sum(Trajectory.length))
    sum_op = sum_op.filter(Trajectory.msm_groups.contains(msmgroup))
    sum_op = sum_op.filter(not_(Trajectory.msm_groups.contains(prev)))
    n_steps = float(sum_op.scalar())
    
    
    p_explore = activation_response(n_new / n_steps, k_factor)

    if len(msmgroup.markov_models) != 2:
        raise ValueError('Only 2 models')
    if not [False, True] == sorted([msm.forcefield.true_kinetics for msm in msmgroup.markov_models]):
        raise ValueError('one needs to be true_kinetics, the other not')
    
    for msm in msmgroup.markov_models:
        if msm.forcefield.true_kinetics:
            msm.model_selection_weight = 1-p_explore
            even_sampling(Session, msm)
        else:
            msm.model_selection_weight = p_explore
            even_sampling(Session, msm)
        
        logger.info("%s selection weight: %f", msm.forcefield.name, msm.model_selection_weight)

#==============================================================================#    


def default(Session, msmgroup):
    """
    If MSMAccelerator can't call any adaptive routine, or none is specified, use
    this one. We default to even sampling.
    """
    even_sampling(Session, msmgroup)


def even_sampling(Session, msmgroup):
    """
    Choose each state with equal probability (uniform dist. over states, ffs)
    """
    for msm in msmgroup.markov_models:
        msm.microstate_selection_weights = np.ones(msmgroup.n_states)
        msm.model_selection_weight = 1


def counts_based(Session, msmgroup, beta=0.0):
    """ 
    This is a simple counts-based adaptive sampling algorithm that places new 
    simulations in states according to how many times those states have been
    "seen".
    
    In [1] it was shown that this algorithm with beta=0 ("mincounts sampling")
    was particularly efficient at exploring space.
    
    Parameters
    ----------
    beta : float (optional)
        A float in [0, inf] that represents an inverse temperature. Low beta
        (< 1.0) corresponds to putting simulations in under-sampled areas, 
        beta = 1.0 corresponds to even sampling, and beta > 1.0 corresponds to
        sampling already sampled regions.
    
    References
    ----------
    ..[1] Weber, J. K. & Pande, V. S. Characterization and Rapid Sampling of 
    Protein Folding Markov State Model Topologies. J. Chem. Theory Comput.
    7, 3405–3411 (2011).
    """

    logger.info("Running counts-based sampling with inverse temp: %f" % beta)

    n_forcefields = len(msmgroup.markov_models)

    microstate_weights = [None] * n_forcefields

    for msm in msmgroup.markov_models:
        sampler = adaptive.CountsSampler(beta)
        msm.microstate_selection_weights = sampler.sample(msm.counts)
        msm.model_selection_weight = 1 # equal chance of selecting any FF
        
        logger.info("%s selection weight: %f", msm.forcefield.name, msm.model_selection_weight)


