import numpy as np
import logging
from msmaccelerator.models import Trajectory, Forcefield, MarkovModel, MSMGroup
from sqlalchemy.sql import and_, or_
logger = logging.getLogger('MSMAccelerator.sampling')

def myfavorite(Session, msmgroup):
    logger.info('Running myfavorite')
    
    for msm in msmgroup.markov_models:
        msm.microstate_selection_weights = np.ones(msmgroup.n_states)
        msm.model_selection_weight = 1
        msm.model_selection_weight
    

#======
def activation_response(x, k):
    """Curve from [0,1] -> [0,1]
    
    """
    if k < 0:
        k = k / (1 - k)
    return x / (1 - k * (x - 1))