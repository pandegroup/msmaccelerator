from __future__ import division
import numpy as np
import logging
from msmaccelerator.models import Trajectory, Forcefield, MarkovModel, MSMGroup
from sqlalchemy.sql import and_, or_
logger = logging.getLogger('MSMAccelerator.sampling')

def myfavorite(Session, msmgroup):
    logger.info('Running myfavorite')
    
    # =============#
    k_factor = 0
    # =============#
           
    prev = Session.query(MSMGroup).get(msmgroup.id - 1)
    if prev is None:
        for msm in msmgroup.markov_models:
            if msm.forcefield.true_kinetics:
                msm.model_selection_weight = 0
            else:
                msm.model_selection_weight = 1
            even_sampling(Session, msm)
        
        
    # number of new states discovered
    n_new = msmgroup.n_states - prev.n_states
    
    # all trajectories in current model
    clause = lambda msm: Trajectory.markov_models.contains(msm)
    t_cur = and_(*[clause(msm) for msm in msmgroup.markov_models])
    # all trajectories in prev model
    t_prev = and_(*[clause(msm) for msm in prev.markov_models])
    
    q = Session.query(Trajectory)
    new_trajectories = q.filter(t_cur).except_(q.filter(t_prev))
    
    # sum of the number of steps in the new trajectories
    steps = sum(t.length - 1 for t in new_trajectories.values(Trajectory.length))
    
    p_exploit = activation_response(n_new / steps, k_factor)

    if len(msmgroup.markov_models) != 2:
        raise ValueError()
    if not [False, True] == sorted([msm.forcefield.true_kinetics for msm in msmgroup.markov_models]):
        raise ValueError()
    
    for msm in msmgroup.markov_models:
        if msm.forcefield.true_kinetics:
            msm.model_selection_weight = p_exploit
            even(Session, msm)
        else:
            msm.model_selection_weight = 1-p_exploit
            even(Session, msm)
    

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