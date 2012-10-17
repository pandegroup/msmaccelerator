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

import os,sys
import numpy as np
from scipy import io
import shutil
from copy import deepcopy
import subprocess
import logging

import msmbuilder.Trajectory
import cPickle as pickle

from collections import defaultdict
from msmaccelerator import Project
from utils import lru_cache
from models import Trajectory, Forcefield, MarkovModel, MSMGroup
from database import Session, with_db_lock
from utils import load_file, save_file

logger = logging.getLogger('MSMAccelerator.Brain')


def multinomial(weights):
    """Choose randomly from a multinomil distribution
    
    The weights are not required to be normalized
    
    Parameters
    ----------
    weights : np.ndarray
        1D array of weights
        
    Returns
    -------
    index : int
        The index of the selected item
    """
    return np.where(np.random.multinomial(n=1, pvals=weights / np.sum(weights)))[0][0]


@with_db_lock
def generate_job():
    """Generate a single job from the most recent MSMs
    
    No parameters -- reads from the database and from the Project file to get
    info.
    
    Returns
    -------
    traj : models.Trajectory
        An unsaved trajectory. Note that we "attach" the conformation that we want
        to start from to the object as traj.init_pdb.
    """
    
    # get the most recent build round
    msmgroup = Session.query(MSMGroup).order_by(MSMGroup.id.desc()).first()
    if msmgroup is None:
        return _generate_equilibration_job()
    
    # select the forcefield to shoot with
    msm_i = multinomial([msm.model_selection_weight for msm in msmgroup.markov_models])
    model_to_shoot_with = msmgroup.markov_models[msm_i]
    forcefield = model_to_shoot_with.forcefield
    
    # select the microstate to pull from
    selected_microstate = multinomial(model_to_shoot_with.microstate_selection_weights)
    
    # select the conformation to start from uniformly over all the confs
    # in that microstate. Note that these confs count come from any of the
    # forcefields. This is why its critical to have the joint clustering
    confs_per_model = []
    for msm in msmgroup.markov_models:
        n_confs = 0
        
        # it is possible that there is no inverse_assignments if this round
        # was done before any trajs were run in this ff
        # its also possible that this state was only populated by one of the ffs,
        # and sense you can't pickle defaultdicts (they have lambda), msm.inv_assignments
        # is a regular dict
        if os.path.exists(msm.inverse_assignments_fn):
            msm.inv_assignmnets = load_file(msm.inverse_assignments_fn)
            if selected_microstate in msm.inv_assignmnets:
                n_confs = len(msm.inv_assignmnets[selected_microstate])

        confs_per_model.append(n_confs)

    model_to_draw_from = msmgroup.markov_models[multinomial(confs_per_model)]
    trjs, frames = model_to_draw_from.inv_assignmnets[selected_microstate]
    r = np.random.randint(len(trjs))
    #import IPython as ip; ip.embed()
    trj_i, frame_i = trjs[r], frames[r]
    
    conf =  model_to_draw_from.trajectories[trj_i].extract_wet_frame(frame_i)
    
    logger.info('Drew conformation from %s, shooting in %s', model_to_draw_from.forcefield.name, forcefield.name)
    if forcefield == model_to_draw_from.forcefield:
        # we're shooting in the same forcefield the conf was generated in
        t = Trajectory(
            forcefield=forcefield,
            name='Started from round {}, microstate {} (ff {}, trj {}, frame {}) -- continuing in same ff'.format(msmgroup.id, selected_microstate, model_to_draw_from.forcefield.name, trj_i, frame_i),
            mode='Production')

    else:
        # we're shooting in a new forcefield
        # SHOULD WE DO EQUILIBRATION INSTEAD OF PRODUCTION?
        # TJL sayz yes, we definitely should
        # ======
        t = Trajectory(
            forcefield=forcefield,
            name='Started from round {}, microstate {} (ff {}, trj {}, frame {}) -- switching forcfields'.format(msmgroup.id,  selected_microstate, model_to_draw_from.forcefield.name, trj_i, frame_i),
            mode='Production')
        
    t.initialized_from_id = model_to_draw_from.trajectories[trj_i].id
    t.initialized_from_frame = int(frame_i)
        
    # this just "attaches" the conf object to t. the Trajectory model doesn't
    # actually have an init_pdb field -- in the submit method of QMaster we'll
    # save this conf to a file on disk
    t.init_pdb = conf
    return t 
                
        
def _generate_equilibration_job():
    """Generate a single equilibration job from the first forcefield
    
    No parameters -- reads from the database and from the Project file to get info.
    
    Returns
    -------
    traj : models.Trajectory
        An unsaved trajectory. Note that we "attach" the conformation that we want
        to start from to the object as traj.init_pdb.
    """
        
    logger.info('Constructing initial equilibration job')
    conf = msmbuilder.Trajectory.load_from_pdb(Project().pdb_topology_file)
        
        
    if Project().starting_confs_lh5 is None:
        # start from pdb_topology_file
        # copy to a new location so that the 'conf' can be deleted without
        # looseing our topology file
        logger.info('Using pdb topolgy to start equlibration run')
        name = 'equilibration, starting from pdb toplogy'
    else:
        num_frames = msmbuilder.Trajectory.load_from_lhdf(Project().starting_confs_lh5, JustInspect=True)[0]
        r = np.random.randint(num_frames)
        xyz = msmbuilder.Trajectory.read_lhdf_frame(Project().starting_confs_lh5, r)
        conf['XYZList'] = np.array([xyz])
        logger.info('Using frame %s of starting_confs_lh5 (%s) to start equilibration run' % (r, Project().starting_confs_lh5))
        name = 'equilibration, starting from frame %s of starting_confs_lh5 (%s)' %  (r, Project().starting_confs_lh5) 
        
    forcefield = Session.query(Forcefield).first()
    
    trj = Trajectory(forcefield=forcefield,
                     name=name, mode='Equilibration')
    trj.init_pdb = conf
    return trj
        
        
