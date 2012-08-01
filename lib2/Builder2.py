import sys, os
import re
import numpy as np
from scipy import io
import logging

# msmbuilder imports
from msmbuilder.MSMLib import GetCountMatrixFromAssignments
from msmbuilder.MSMLib import EstimateReversibleCountMatrix
from msmbuilder.MSMLib import ErgodicTrim, ApplyMappingToAssignments
from msmbuilder import clustering, metrics
import msmbuilder.Trajectory
import msmbuilder.Serializer
from msmbuilder.assigning import assign_in_memory

from models import Trajectory, Forcefield, MarkovModel, MSMGroup
from database import Session
logger = logging.getLogger('MSMAccelerator.Builder')

class Builder(object):
    """
    The Builder class is responsable for managing the construction of MSMs.
    it takes in the data generated by the QMaster and adds it to the database
    """
    
    def __init__(self, project):
        self.project = project
    
    @property
    def n_rounds(self):
        return Session.query(MSMGroup).count()
        
    def is_sufficient_new_data(self):
        """Is there sufficient new data to build a new round?
        
        Returns
        -------
        truth : boolean
            True if there is sufficient new data for a new round
        """
        
        # get most recent set of MSMs built
        msmgroup = Session.query(MSMGroup).order_by(MSMGroup.id.desc()).first()
        if msmgroup is None:
            n_built = 0
        else:        
            # the number of unique trajectories that are part of this msmgroup
            # by constructing a query for the union of trajectories in any of the msms
            # in the msm group
            q = Session.query(Trajectory)
            for msm in msmgroup.markov_models:
                q2 = Session.query(Trajectory)
                q2.filter(Trajectory.markov_models.contains(msm))
                q2.filter(Trajectory.returned_time != None)
                q.union(q2)
            n_built = q.count()
            print q.all()
            ei
        
        # number of trajs in the database
        n_total = Session.query(Trajectory).filter(Trajectory.returned_time != None).count()
        
        truth = n_total > n_built + self.project.num_trajs_sufficient_for_round
        
        logger.info("%d trajs total, %d trajs built. Sufficient? %s", n_total, n_built, truth)
        return truth
    
    
    def run_round(self, checkdata=True):
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
            if not self.is_sufficient_new_data():
                return False
        else:
            logger.info("Skipping check for adequate data.")
        
        # use all the data together to get the cluster centers
        generators = self.joint_clustering()
        
        msmgroup = MSMGroup()
        for ff in Session.query(Forcefield).all():
            msm = self.build_msm(ff, generators=generators)
            msmgroup.markov_models.append(msm)
        
        Session.add(msmgroup)
        Session.flush()
        logger.info("Round completed sucessfully")
        
        return True
        
        
    def joint_clustering(self):
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
        load = lambda v: msmbuilder.Trajectory.LoadTrajectoryFile(v)
        loaded_trjs = [load(t.lh5_fn)[::self.project.stride] for t in db_trajs]
        
        clusterer = self.project.clusterer(trajectories=loaded_trjs)
        return clusterer.get_generators_as_traj()

            
    def build_msm(self, forcefield, generators):
        """Build an MSM for this forcefield using the most recent trajectories
        in the database.
        
        If supplied, use the supplied generators
        
        Parameters
        ----------
        forcefield : models.Forcefield
            database entry on the forcefield that we build for
        generators : msmbuilder.Trajectory

        Returns
        -------
        msm : models.MarkovModel
        
        """
        
        # I want to use assign_in_memory, which requires an msmbuilder.Project
        # so, lets spoof it
        trajs = Session.query(Trajectory).filter(Trajectory.returned_time != None).all()
        
        class BuilderProject(object):
            def __init__(self):
                self['NumTrajs'] = len(trajs)
                self['TrajLengths'] = np.array([t.length for t in trajs])
            
            def LoadTraj(self, trj_index):
                if trj_index < 0 or trj_index > len(trajs):
                    raise IndexError('Sorry')
                return msmbuilder.Trajectory.LoadTrajectoryFile(trajs[trj_index].lh5_fn)
        
        
        logger.info('Assigning...')
        assignments, distances = assign_in_memory(self.project.metric, generators, BuilderProject())
        
        logger.info('Getting counts...')
        counts = self.construct_counts_matrix(assignments)
        
        return MarkovModel(counts=counts, assignments=assignments, distances=distances,
            forcefield=forcefield, trajectories=trajs)
            
        
        
    def construct_counts_matrix(self, assignments):
        """Build and return a counts matrix from assignments. Symmetrize either
        with transpose or MLE based on the value of the self.symmetrize variable
        
        Also modifies the assignments file that you pass it to reflect ergodic
        trimming"""
        
        n_states  = np.max(assignments.flatten()) + 1
        raw_counts = GetCountMatrixFromAssignments(assignments, n_states,
                                    LagTime=self.project.lagtime, Slide=True)
        
        ergodic_counts = None
        if self.project.trim:
            raise NotImplementedError(('Trimming is not yet supported because '
                'we need to keep track of the mapping from trimmed to '
                ' untrimmed states for joint clustering to be right'))
            try:
                ergodic_counts, mapping = ErgodicTrim(raw_counts)
                ApplyMappingToAssignments(assignments, mapping)
                counts = ergodic_counts
            except Exception as e:
                logger.warning("ErgodicTrim failed with message '{0}'".format(e))
        else:
            logger.info("Ignoring ergodic trimming")
            counts = raw_counts
        
        if self.project.symmetrize == 'transpose':
            logger.debug('Transpose symmetrizing')
            counts = counts + counts.T
        elif self.project.symmetrize == 'mle':
            logger.debug('MLE symmetrizing')
            counts = EstimateReversibleCountMatrix(counts)
        elif self.project.symmetrize == 'none' or self.project.symmetrize == None:
            logger.debug('Skipping symmetrization')
        else:
            raise ValueError("Could not understand symmetrization method: %s" % self.project.symmetrize)
        
        return counts

