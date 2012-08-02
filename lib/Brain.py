import os,sys
import numpy as np
from scipy import io
import shutil
from copy import deepcopy
import subprocess
import logging

import msmbuilder.Trajectory
import msmbuilder.Serializer
import cPickle as pickle

from collections import defaultdict
import Sampling
from utils import lru_cache
from models import Trajectory, Forcefield, MarkovModel, MSMGroup
from database import Session
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


class Brain(object):
    def __init__(self, project):
        self.project = project
    
    def generate_job(self):
        """Generate a single job from the most recent MSMs
        
        Parameters
        ----------
        
        Returns
        -------
        An unsaved models.Job db object
        """
        
        # get the most recent build round
        msmgroup = Session.query(MSMGroup).order_by(MSMGroup.id.desc()).first()
        if msmgroup is None:
            return self._generate_equilibration_job()
        
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
            with open(msm.inverse_assignments_fn) as f:
                msm.inv_assignmnets = pickle.load(f)
            n_confs = len(msm.inv_assignmnets[selected_microstate])
            confs_per_model.append(n_confs)

        model_to_draw_from = msmgroup.markov_models[multinomial(confs_per_model)]
        assigned_to = model_to_draw_from.inv_assignmnets[selected_microstate]
        r = np.random.randint(len(assigned_to))
        trj_i, frame_i = assigned_to[r]
        
        conf =  model_to_draw_from.trajectories[trj_i].extract_wet_frame(frame_i)
        
        
        logger.info('Drew conformation from %s, shooting in %s', model_to_draw_from.forcefield.name, forcefield.name)
        if forcefield == model_to_draw_from.forcefield:
            # we're shooting in the same forcefield the conf was generated in
            t = Trajectory(
                forcefield=forcefield,
                name='Started from round {}, microstate {} (ff {}, trj {}, frame {}) -- continuing in same ff'.format(msmgroup.id,
                    selected_microstate, model_to_draw_from.forcefield.name, trj_i, frame_i),
                mode='Production')
        else:
            # we're shooting in a new forcefield
            # NOTE
            # ======
            # SHOULD WE DO EQUILIBRATION INSTEAD OF PRODUCTION?
            # ======
            t = Trajectory(
                forcefield=forcefield,
                name='Started from round {}, microstate {} (ff {}, trj {}, frame {}) -- switching forcfields'.format(msmgroup.id,
                    selected_microstate, model_to_draw_from.forcefield.name, trj_i, frame_i),
                mode='Production')
        
        t.initialized_from_id = model_to_draw_from.trajectories[trj_i].id
        t.initialized_from_frame = int(frame_i)
        
        # this just "attaches" the conf object to t. the Trajectory model doesn't
        # actually have an init_pdb field -- in the submit method of QMaster we'll
        # save this conf to a file on disk
        t.init_pdb = conf
        return t 
                
        
    def _generate_equilibration_job(self):
        """Generate a single equilibration job from the first forcefield
        
        Parameters
        ----------
        tmp_dir : str
            Location to save PDBs that will be used as initial structures for
            the job
        
        Returns
        -------
        trajectory : models.Trajectory
        
        """
        
        
        logger.info('Constructing initial equilibration job')
        conf = msmbuilder.Trajectory.LoadFromPDB(self.project.pdb_topology_file)
        
        
        if self.project.starting_confs_lh5 is None:
            # start from pdb_topology_file
            # copy to a new location so that the 'conf' can be deleted without
            # looseing our topology file
            logger.info('Using pdb topolgy to start equlibration run')
            name = 'equilibration, starting from pdb toplogy'
        else:
            num_frames = msmbuilder.Trajectory.LoadFromLHDF(self.project.starting_confs_lh5, JustInspect=True)[0]
            r = np.random.randint(num_frames)
            xyz = msmbuilder.Trajectory.ReadLHDF5Frame(self.project.starting_confs_lh5, r)
            conf['XYZList'] = np.array([xyz])
            logger.info('Using frame %s of starting_confs_lh5 (%s) to start equilibration run' % (r, self.project.starting_confs_lh5))
            name = 'equilibration, starting from frame %s of starting_confs_lh5 (%s)' %  (r, self.project.starting_confs_lh5) 
        
        forcefield = Session.query(Forcefield).first()
        
        trj = Trajectory(forcefield=forcefield,
            name=name, mode='Equilibration')
        trj.init_pdb = conf
        return trj
        
        