import IPython as ip
import tempfile
import scipy.io
import scipy.sparse
import numpy as np

from models import Trajectory, Forcefield, MarkovModel, MSMGroup
from msmaccelerator import Project
project = Project('rmcgibbo/project.yaml')
db = project.connect_to_db()


def multinomial(weights):
    return np.where(np.random.multinomial(n=1, pvals=weights))[0][0]
    
def setup():
    t1 = Trajectory(name='t1')
    t2 = Trajectory(name='t2')
    t3 = Trajectory(name='t3')

        
    group = MSMGroup()
    amber = Forcefield(name='amber')
    msm1 = MarkovModel(msm_group=group, trajectories=[t1,], forcefield=amber)
    msm1.model_selection_weight = 0.5

    #msm1.inverse_assignments = {'a': 1}
    msm2 = MarkovModel(msm_group=group, trajectories=[t2,], forcefield=amber)
    msm2.model_selection_weight = 0.5
    #msm2.inverse_assignments =  {'a': 2}
    group.microstate_selection_weights = np.array([1,2,3,4,5])
    
    
    db.add_all([t1, t2, t3, group, msm1, msm2])


if __name__ == '__main__':    
    setup()
    m = db.query(MarkovModel).first()
    ip.embed()

