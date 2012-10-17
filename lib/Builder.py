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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys, os
import numpy as np
import logging
import numbers
from collections import defaultdict

# msmbuilder imports
from msmbuilder import MSMLib
from msmbuilder import clustering
from msmbuilder import metrics
import msmbuilder.Trajectory
from msmbuilder.assigning import assign_in_memory

from sqlalchemy.sql import and_, or_
from msmaccelerator import Project
from models import Trajectory, Forcefield, MarkovModel, MSMGroup
from database import Session, with_db_lock
import sampling
from utils import load_file, save_file

logger = logging.getLogger('MSMAccelerator.Builder')

@with_db_lock
def n_rounds():
    "Number of groups of MSMs that have been built"
    return Session.query(MSMGroup).count()

@with_db_lock
def is_sufficient_new_data():
    """Is there sufficient new data to build a new round?
        
    Returns
    -------
    truth : boolean
        True if there is sufficient new data for a new round
    """

    qg = Session.query(MSMGroup)
    qt = Session.query(Trajectory)
        
    msmgroup = qg.order_by(MSMGroup.id.desc()).first()
    if msmgroup is not None:
        n_built = qt.filter(Trajectory.msm_groups.contains(msmgroup)).count()
    else:
        n_built = 0
        
    n_total = qt.filter(Trajectory.returned_time != None).count()
        
    truth = n_total >= n_built + Project().num_trajs_sufficient_for_round
        
    logger.info("%d trajs total, %d trajs built. Sufficient? %s", n_total, n_built, truth)
    return truth


@with_db_lock
def run_round(checkdata=True):
    """Activate the builder and build new MSMs (if necessary)
    
    First, check to see if there is enough data are to warrant building a
    new set of MSMs. Assuming yes, do a joint clustering over all of the
    data, and then build MSMs for each forcefield on that state space.
    
    Parameters
    ----------
    checkdata : boolean, optional
         If False, skip the checking process
    
    Returns
    -------
    happened : boolean
        True if we actually did a round of MSM building, False otherwise
    """
        
    if checkdata:
        logger.info("Checking if sufficient data has been acquired.")
        if not is_sufficient_new_data():
            return False
    else:
        logger.info("Skipping check for adequate data.")
        
    # use all the data together to get the cluster centers
    generators, db_trajs = joint_clustering()
        
    msmgroup = MSMGroup(trajectories=db_trajs)
    for ff in Session.query(Forcefield).all():
        trajs = filter(lambda t: t.forcefield == ff, db_trajs)
        msm = build_msm(ff, generators=generators, trajs=trajs)
        msmgroup.markov_models.append(msm)
        
    # add generators to msmgroup
    Session.add(msmgroup)
    Session.flush()
    msmgroup.populate_default_filenames()
        
    msmgroup.trajectories = db_trajs
    msmgroup.n_states = len(generators)
    save_file(msmgroup.generators_fn, generators)

    for msm in msmgroup.markov_models:
        msm.populate_default_filenames()
        if hasattr(msm, 'counts'):
            save_file(msm.counts_fn, msm.counts)
        if hasattr(msm, 'assignments'):
            save_file(msm.assignments_fn, msm.assignments)
        if hasattr(msm, 'distances'):
            save_file(msm.distances_fn, msm.distances)
            save_file(msm.inverse_assignments_fn, dict(MSMLib.invert_assignments(msm.assignments)))
        
    # ====================================================================== #
    # HERE IS WHERE THE ADAPTIVE SAMPLING ALGORITHMS GET CALLED
    # The obligation of the adaptive_sampling routine is to set the 
    # model_selection_weight on each MSM/forcefield and the microstate selection
    # weights
    
    try:
        Project().adaptive_sampling(Session, msmgroup, **Project().adaptive_parameters)
        
        for msm in msmgroup.markov_models:
            if not isinstance(msm.model_selection_weight, numbers.Number):
                raise ValueError('model selection weight on %s not set correctly' % msm)
            if not isinstance(msm.microstate_selection_weights, np.ndarray):
                raise ValueError('microstate_selection_weights on %s not set correctly' % msm)
                
    except Exception as e:
        logging.error('ADAPTIVE SAMPLING ERROR -- Could not run sampling algorithm')
        logging.error(e)
        logging.error('APPLYING DEFAULT ADAPTIVE ALGORITHM')
        sampling.default(Session, msmgroup)
        
    #======================================================================= #

        
    Session.flush()       
    logger.info("Round completed sucessfully")
    return True

@with_db_lock
def joint_clustering():
    """Jointly cluster the the data from all of the forcefields
    
    Returns
    -------
    generators : msmbuilder.Trajectory
    """
    logger.info('Running joint clustering')
    
    # load up all the trajs in the database
    db_trajs = Session.query(Trajectory).filter(Trajectory.returned_time != None).all()
    if len(db_trajs) == 0:
        raise RuntimeError()
        
    # load the xyz coordinates from disk for each trajectory
    load = lambda v: msmbuilder.Trajectory.load_trajectory_file(v)
    loaded_trjs = [load(t.lh5_fn)[::Project().stride] for t in db_trajs]
    
    clusterer = Project().clusterer(trajectories=loaded_trjs)
    return clusterer.get_generators_as_traj(), db_trajs

            
def build_msm(forcefield, generators, trajs):
    """Build an MSM for this forcefield using the most recent trajectories
    in the database.
        
    If supplied, use the supplied generators
        
    Parameters
    ----------
    forcefield : models.Forcefield
        database entry on the forcefield that we build for
    generators : msmbuilder.Trajectory
    trajs : list of models.Trajectory

    Returns
    -------
    msm : models.MarkovModel
    
    """
        
    # I want to use assign_in_memory, which requires an msmbuilder.Project
    # so, lets spoof it
        
        
    if len(trajs) == 0:
        return MarkovModel(forcefield=forcefield)
        
    class BuilderProject(dict):
        @property
        def n_trajs(self):
            return len(trajs)

        @property
        def traj_lengths(self):
            return np.array([t.length for t in trajs])

        def load_traj(self, trj_index):
            if trj_index < 0 or trj_index > len(trajs):
                raise IndexError('Sorry')
            return msmbuilder.Trajectory.load_trajectory_file(trajs[trj_index].lh5_fn)
        
        
    logger.info('Assigning...')
    assignments, distances = assign_in_memory(Project().metric, generators, BuilderProject())
        
    logger.info('Getting counts...')
    counts = construct_counts_matrix(assignments)
        
        
    model = MarkovModel(forcefield=forcefield, trajectories=trajs)
    model.counts = counts
    model.assignments = assignments
    model.distances = distances
        
    return model
        
        
def construct_counts_matrix(assignments):
    """Build and return a counts matrix from assignments.
    
    Symmetrize either with transpose or MLE based on the value of the
    self.symmetrize variable
        
    Also modifies the assignments file that you pass it to reflect ergodic
    trimming
    
    Parameters
    ----------
    assignments : np.ndarray
        2D array of MSMBuilder assignments
    
    Returns
    -------
    counts : scipy.sparse.csr_matrix
        transition counts
    
    """
        
    n_states  = np.max(assignments.flatten()) + 1
    raw_counts = MSMLib.get_count_matrix_from_assignments(assignments, n_states,
                                               lag_time=Project().lagtime,
                                               sliding_window=True)
        
    ergodic_counts = None
    if Project().trim:
        raise NotImplementedError(('Trimming is not yet supported because '
                                   'we need to keep track of the mapping from trimmed to '
                                   ' untrimmed states for joint clustering to be right'))
        try:
            ergodic_counts, mapping = MSMLib.ergodic_trim(raw_counts)
            MSMLib.apply_mapping_to_assignments(assignments, mapping)
            counts = ergodic_counts
        except Exception as e:
            logger.warning("MSMLib.ergodic_trim failed with message '{0}'".format(e))

    else:
        logger.info("Ignoring ergodic trimming")
        counts = raw_counts
        
    if Project().symmetrize == 'transpose':
        logger.debug('Transpose symmetrizing')
        counts = counts + counts.T
    elif Project().symmetrize == 'mle':
        logger.debug('MLE symmetrizing')
        counts = MSMLib.mle_reversible_count_matrix(counts)
    elif Project().symmetrize == 'none' or (not Project().symmetrize):
        logger.debug('Skipping symmetrization')
    else:
        raise ValueError("Could not understand symmetrization method: %s" % Project().symmetrize)
        
    return counts
