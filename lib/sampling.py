import numpy as np
import logging
from models import Trajectory, Forcefield, MarkovModel, MSMGroup
logger = logging.getLogger('MSMAccelerator.sampling')

def myfavorite(Session, msmgroup):
    logger.info('Running myfavorite')
    
    # get number of states in previous rounds' models
    prev = Session.query(MSMGroup).get(msmgroup.id - 1):
    if prev is not None:
        n_prev_states = prev.n_states
    else:
        n_prev_state = 0
    
    # fraction of the current states which were in previous model
    # f= 0 if all the states are new -- not converged
    # f= 1 if we're converged
    f = n_prev_states / msmgroup.n_states
    
    
    