import os,sys
import numpy as np
from scipy import io
import shutil
import Sampling
from msmbuilder.Trajectory import Trajectory
from msmbuilder import Serializer
from collections import defaultdict
from utils import lru_cache
from copy import deepcopy
import subprocess
#import IPython as ip
import logging

def multinomial(weights):
    return np.where(np.random.multinomial(n=1, pvals=weights))[0][0]

class Brain(object):
    def __init__(self, project):
        self.logger = logging.getLogger('MSMAccelerator.Brain')
        self.project = project
        self.round_aux_data = {}
        
        # overwrite method
        self.adaptive_sampling = getattr(Sampling, project.adaptive_sampling)
    
    @staticmethod
    def adaptive_sampling(ff_data, num_jobs):
        raise NotImplementedError(('You need to provide an adaptive sampling algorithm '
            'which will be injected into this class at instantiation'))

    def generate_job(self, round_num, tmp_dir='/tmp'):
        """Generate a single job from the MSM built in round_num"""
        
        if round_num < 0:
            return self._generate_initial_job(tmp_dir=tmp_dir)
        
        ffs, ff_weights, microstate_weights = self._aux_data(round_num)
        
        # choose which forcefield based on a multinomial
        ff_index = multinomial(ff_weights)
        # choose which microstate within that forcefield
        microstate = multinomial(microstate_weights[ff_index])
        
        ff = ffs[ff_index]
        assigned_to = ff['inverse_assignments'][microstate]
        r = np.random.randint(len(assigned_to))
        trj_i, frame_i = assigned_to[r]
        conf = self.extract_xtc_frame(ff, trj_i, frame_i, tmp_dir=tmp_dir)
        return {'conf': conf,
                'name': 'Started from round {0}, microstate {1} (trj {2}, frame {3})'.format(round_num, microstate, trj_i, frame_i),
                'ff': ff['name'],
                'water': ff['water'],
                'driver': os.path.abspath(ff['driver']),
                'threads': ff['threads'],
                'mode': 'Production',
                'output_extension': ff['output_extension']}
    
    @lru_cache(1)
    def _aux_data(self, round_num):
        """Loads up the inverse assignments, counts matricies and adaptive sampling
        multinomial weights for the given round. Since this function is memoized, it
        doesn't have to be called very often"""
        
        self.logger.info('Loading assignments, counts for round {0} into memory'.format(round_num))
        
        loaded_ffs = deepcopy(self.project.forcefields)
        for forcefield in loaded_ffs:
            ff_name = forcefield['name']

            assignments = Serializer.LoadData(self.project.assignments_fn(ff_name, round_num))
            inverse_assignments = defaultdict(lambda:  [])
            for i in xrange(assignments.shape[0]):
                for j in xrange(assignments.shape[1]):
                    if assignments[i,j] != -1:
                        inverse_assignments[assignments[i,j]].append((i,j))
                        
            forcefield['counts'] = io.mmread(self.project.counts_fn(ff_name, round_num))
            forcefield['inverse_assignments'] = inverse_assignments
        
        ff_weights, microstate_weights = self.adaptive_sampling(loaded_ffs)

        return loaded_ffs, ff_weights, microstate_weights
                
    @deprecated
    def extract_xtc_frame(self, forcefield, trj_num, frame_index, tmp_dir='/tmp'):
        """Extract a PDB of the solvated xtc from a certain forcefield, traj_num, frame_index pair
        Returns the path to the PDB"""
        
        tt = self.project.traj_log
        traj = [t for t in tt if t['lh5_trajnum'] == trj_num and t['ff'] == forcefield['name']][0]
        
        if self.project.method == 'explicit':
            xtc_path = traj['wet_xtc']
            conf = Trajectory.LoadTrajectoryFile(traj['last_wet_snapshot'])
        else:
            conf = Trajectory.LoadTrajectoryFile(self.project.pdb_topology_file)
            xtc_path = traj['dry_xtc']
        
        xyz = Trajectory.ReadFrame(xtc_path, frame_index)
        if not (xyz.shape == conf['XYZList'][0,:,:].shape):
            raise Exception("Number of atoms is wrong: xyz.shape={0}, conf['XYZList'][0,:,:].shape={1}".format(xyz.shape, conf['XYZList'][0,:,:].shape))
           
        conf['XYZList'] = np.array([xyz])

        output = '%s/%d.pdb' % (tmp_dir, np.random.randint(sys.maxint))
        conf.SaveToPDB(output)
        
        self.logger.info('Ripping out frame {frame_ind} from {xtc} trjconv'.format(frame_ind=frame_index, xtc=xtc_path))
        return output
            
    def _generate_initial_job(self, tmp_dir='/tmp'):
        self.logger.info('Constructing equlilibration job')
        ff = self.project.forcefields[0]
        
        conf_fn = '%s/%d.pdb' % (tmp_dir, np.random.randint(sys.maxint))
        
        if self.project.starting_confs_lh5 is None:
            # start from pdb_topology_file
            # copy to a new location so that the 'conf' can be deleted without
            # looseing our topology file
            shutil.copy(self.project.pdb_topology_file, conf_fn)
            self.logger.info('Using pdb topolgy to start equlibration run')
            name = 'equilibration, starting from pdb toplogy'
        else:
            conf = Trajectory.LoadFromPDB(self.project.pdb_topology_file)
            num_frames = Trajectory.LoadFromLHDF(self.project.starting_confs_lh5, JustInspect=True)[0]
            r = np.random.randint(num_frames)
            xyz = Trajectory.ReadLHDF5Frame(self.project.starting_confs_lh5, r)
            conf['XYZList'] = np.array([xyz])
            conf.SaveToPDB(conf_fn)
            self.logger.info('Using frame %s of starting_confs_lh5 (%s) to start equilibration run' % (r, self.project.starting_confs_lh5))
            name = 'equilibration, starting from frame %s of starting_confs_lh5 (%s)' %  (r, self.project.starting_confs_lh5) 
            
        return {
            'name': name,
            'ff': ff['name'],
            'water': ff['water'],
            'driver': os.path.abspath(ff['driver']),
            'threads': ff['threads'],
            'mode': 'Equilibration',
            'conf': conf_fn,
            'output_extension': ff['output_extension']
            }
