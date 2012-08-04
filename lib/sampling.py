from __future__ import division
import numpy as np
import logging
import IPython as ip
from msmaccelerator.models import Trajectory, Forcefield, MarkovModel, MSMGroup
from sqlalchemy.sql import and_, or_, not_
from sqlalchemy import func


logger = logging.getLogger('MSMAccelerator.sampling')

def myfavorite(Session, msmgroup):
    logger.info('Running myfavorite')
    
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
    

def default(Session, msmgroup):
    for msm in msmgroup.markov_models:
        msm.microstate_selection_weights = np.ones(msmgroup.n_states)
        msm.model_selection_weight = 1


#==============================================================================#

def even_sampling(Session, msm):
    msm.microstate_selection_weights = np.ones(msm.msm_group.n_states)


    #==============================================================================#

def activation_response(x, k):
    """Curve from [0,1] -> [0,1]
    
    """
    if k < 0:
        k = k / (1 - k)
    return x / (1 - k * (x - 1))
