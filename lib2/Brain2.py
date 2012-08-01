import os,sys
import numpy as np
from scipy import io
import shutil
from copy import deepcopy
import subprocess
import logging

import msmbuilder.Trajectory
import msmbuilder.Serializer

from collections import defaultdict
import Sampling
from utils import lru_cache
from models import Trajectory, Forcefield, MarkovModel, MSMGroup

def multinomial(weights):
    return np.where(np.random.multinomial(n=1, pvals=weights))[0][0]


class Brain(object):
    def __init__(self, project):
        self.logger = logging.getLogger('MSMAccelerator.Brain')
        self.project = project
        self.db = project.db
        
    #     # overwrite method
    #     self.adaptive_sampling = getattr(Sampling, project.adaptive_sampling)
    # 
    # @staticmethod
    # def adaptive_sampling(ff_data, num_jobs):
    #     raise NotImplementedError(('You need to provide an adaptive sampling algorithm '
    #         'which will be injected into this class at instantiation'))
    
    
    def generate_job(self):
        """Generate a single job from the most recent MSMs
        
        Parameters
        ----------
        
        Returns
        -------
        An unsaved models.Job db object
        """
        
        # get the most recent build round
        msmgroup = self.db.query(MSMGroup).order_by(MSMGroup.id.desc()).first()
        if msmgroup is None:
            return self._generate_equilibration_job()
        
        # select the microstate from multinomial associated with the group
        # of models, all on the same statespace
        selected_microstate = multinomial(msmgroup.microstate_selection_weights)
        
        # select the conformation to start from uniformly over all the confs
        # in that microstate. Note that these confs count come from any of the
        # forcefields
        confs_per_model = []
        for msm in msm_group.markov_models:
            confs_per_model.append(len(msm.inverse_assignmnents[selected_microstate]))
        model_to_draw_from = msmgroup.markov_models[multinomial(confs_per_model)]
        assigned_to = model_to_draw_from.inverse_assignments[selected_microstate]
        r = np.random.randint(len(assigned_to))
        trj_i, frame_i = assigned_to[r]
        
        conf =  model_to_draw_from.trajectories[trj_i].extract_wet_frame(frame_i)
        
        # select the forcefield to shoot with
        msm_i = multinomial([msm.model_selection_weight for msm in msm_group.markov_models])
        forcefield = msmgroup.markov_models[msm_i].forcefield
            
        if forcefield == model_to_draw_from.forcefield:
            # we're shooting in the same forcefield the conf was generated in
            return Trajectory(
                init_pdb=conf,
                forcefield=forcefield,
                name='Started from round {}, microstate {} (ff {}, trj {}, frame {}) -- continuing in same ff'.format(build_round.id,
                    selected_microstate, model_to_draw_from.forcefield.name, trj_i, frame_i),
                mode='Production')
        else:
            # we're shooting in a new forcefield
            # NOTE
            # ======
            # SHOULD WE DO EQUILIBRATION INSTEAD OF PRODUCTION?
            # ======
            return Trajectory(
                init_pdb=conf,
                forcefield=forcefield,
                name='Started from round {}, microstate {} (ff {}, trj {}, frame {}) -- switching forcfields'.format(build_round.id,
                    selected_microstate, model_to_draw_from.forcefield.name, trj_i, frame_i),
                mode='Production')
                
        
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
        
        
        self.logger.info('Constructing initial job')
        conf = msmbuilder.Trajectory.LoadFromPDB(self.project.pdb_topology_file)
        
        
        if self.project.starting_confs_lh5 is None:
            # start from pdb_topology_file
            # copy to a new location so that the 'conf' can be deleted without
            # looseing our topology file
            self.logger.info('Using pdb topolgy to start equlibration run')
            name = 'equilibration, starting from pdb toplogy'
        else:
            num_frames = msmbuilder.Trajectory.LoadFromLHDF(self.project.starting_confs_lh5, JustInspect=True)[0]
            r = np.random.randint(num_frames)
            xyz = msmbuilder.Trajectory.ReadLHDF5Frame(self.project.starting_confs_lh5, r)
            conf['XYZList'] = np.array([xyz])
            self.logger.info('Using frame %s of starting_confs_lh5 (%s) to start equilibration run' % (r, self.project.starting_confs_lh5))
            name = 'equilibration, starting from frame %s of starting_confs_lh5 (%s)' %  (r, self.project.starting_confs_lh5) 
        
        forcefield = self.db.query(Forcefield).first()
        
        trj = Trajectory(forcefield=forcefield,
            name=name, mode='Equilibration')
        trj.init_pdb = conf
        
        return trj
        
        