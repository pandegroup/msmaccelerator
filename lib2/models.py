# DECLARE MODELS
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import (Column, Integer, String, DateTime,
                        Float, ForeignKey, Table, PickleType)
Base = declarative_base()
STRING_LEN = 500

import msmbuilder.Trajectory
from database import file_store_metaclass
from Project2 import Project

class Forcefield(Base):
    """
    database model describing a forcefield and how to run simulations using
    it (driver, threads, etc)
    """
    __tablename__ = 'forcefields'
    
    id = Column(Integer, primary_key=True)
    name = Column(String(STRING_LEN), unique=True)
    water = Column(String(STRING_LEN))
    driver = Column(String(STRING_LEN))
    threads = Column(Integer)
    output_extension = Column(String(STRING_LEN))
        
    def __repr__(self):
        return "<Forcefield(name={}, water={}, driver={})>".format(self.name,
            self.water, self.driver)


from sqlalchemy.orm.interfaces import MapperExtension
from sqlalchemy.orm.attributes import InstrumentedAttribute
class MAPEX(MapperExtension):
    @staticmethod
    def after_insert(mapper, connection, instance):
        for k,v in instance.__class__.__dict__.items():
            if k.startswith('_') and k.endswith('_fn') and isinstance(v, InstrumentedAttribute):
                if getattr(instance, k) is None:
                    name_of_default_method = '_default' + k + '_'
                    value = getattr(instance, name_of_default_method)()
                    setattr(instance, k, value)

class Trajectory(Base):
    """
    database model describing a trajectory
    """
    
    __tablename__ = 'trajectories'
    __metaclass__ = file_store_metaclass
    __mapper_args__ = {'extension' : MAPEX}
    
    id = Column(Integer, primary_key=True)
    _init_pdb_fn = Column(String(STRING_LEN))
    _wqlog_fn = Column(String(STRING_LEN))
    _dry_xtc_fn = Column(String(STRING_LEN))
    _wet_xtc_fn = Column(String(STRING_LEN))
    _lh5_fn = Column(String(STRING_LEN))
    _last_wet_snapshot_fn = Column(String(STRING_LEN))

    host = Column(String(STRING_LEN))
    mode = Column(String(STRING_LEN))
    name = Column(String(STRING_LEN))
    returned_time = Column(Float)
    submit_time = Column(Float)
    traj_length = Column(Integer) # in frames
    
    # many to one relationship
    forcefield_id = Column(Integer, ForeignKey('forcefields.id'))
    forcefield = relationship("Forcefield", backref=backref('trajectories', order_by=id))
    
    
    def __basefn(self):
            return '{base}/{ff}/'.format(base=Project.instance.project_dir, ff=self.forcefield.name)
    def _default_init_pdb_fn_(self):
        return self.__basefn() + 'Inits/{id}.pdb'.format(id=self.id)
    def _default_wqlog_fn_(self):
        return self.__basefn() + 'WQLogs/{id}.log.txt'.format(id=self.id)
    def _default_dry_xtc_fn_(self):
        return self.__basefn() + 'DryXTCs/{id}.xtc'.format(id=self.id)
    def _default_wet_xtc_fn_(self):
        return self.__basefn() + 'WetXTCs/{id}.xtc'.format(id=self.id)
    def _default_lh5_fn_(self):
        return self.__basefn() + 'LH5s/{id}.lh5'.format(id=self.id)
    def _default_last_wet_snapshot_fn_(self):
        return self.__basefn() + 'LastWetSnapshots/{id}.lh5'.format(id=self.id)
    




    def __repr__(self):
        return "<Trajectory(name={})>".format(self.name)
        
    def extract_wet_frame(self, frame_index):
        """Extract a PDB of a solvated structure (single frame) in this trajectory
        
        Parameters
        ----------
        frame_index : int
            Index of the frame to extract
        
        Returns
        -------
        path : msmbuilder.Trajectory
            
        """
        
        if self.wet_xtc_fn is not None:
            conf = msmbuilder.Trajectory.LoadTrajectoryFile(self.last_wet_snapshot)
            xyz = msmbuilder.Trajectory.ReadFrame(self.wet_xtc_fn, frame_index)
        else:
            conf = msmbuilder.Trajectory.LoadTrajectoryFile(self.project.pdb_topology_file)
            xyz = msmbuilder.Trajectory.ReadFrame(xtc_path, frame_index)

        if not (xyz.shape == conf['XYZList'][0,:,:].shape):
            raise Exception("Number of atoms is wrong: xyz.shape={0}, conf['XYZList'][0,:,:].shape={1}".format(xyz.shape, conf['XYZList'][0,:,:].shape))

        conf['XYZList'] = np.array([xyz])

        return conf



class MarkovModel(Base):
    """
    database model describing an MSM built from a set of trajectories, all in the
    same forcefield
    
    Each BuildRound contains multiple MarkovModels
    """
    __tablename__ = 'markov_models'
    __metaclass__ = file_store_metaclass
    
    id = Column(Integer, primary_key=True)
    _counts_fn = Column(String(STRING_LEN))
    _assignments_fn = Column(String(STRING_LEN))
    _inverse_assignments_fn = Column(String(STRING_LEN))
    
    def __basefn(self):
         return '{base}/{ff}/models/{id}/'.format(base=Project.instance.project_dir,
            ff=self.forcefield.name, id=self.msm_group_id)
    def _default_counts_fn_(self):
        return self.__basefn() + 'tCounts.mtx'
    def _default_inverse_assignments_fn_(self):
        return self.__basefn() + 'inverse_assignments.pickl'
    def _default_assignments_fn_(self):
        return self.__basefn() + 'assignments.h5'
    
    forcefield_id = Column(Integer, ForeignKey('forcefields.id'))
    forcefield = relationship("Forcefield", backref=backref('markov_models',
        order_by=id))
    
    # a many to many relationship between models and trajectories -- each model
    # can be have multiple trajectories and each trajectory can be in multiple
    # models
    trajectories = relationship("Trajectory",  backref='markov_models',
        secondary=Table('markov_model_trajectories', Base.metadata,
                Column('markov_model_id', Integer, ForeignKey('markov_models.id')),
                Column('trajectory_id', Integer, ForeignKey('trajectories.id'))))
    
    msm_group_id = Column(Integer, ForeignKey('msm_groups.id'))
    msm_group = relationship('MSMGroup', backref=backref('markov_models',
        order_by=id))
    
    # When selecting which forcefield to shoot from after a round, choose based
    # on a multinomial with this weight
    model_selection_weight = Column(Float)
    

    
    def __repr__(self):
        return "<MarkovModel(id={}, msm_group={}, forcefield={})>".format(self.id,
            self.msm_group_id, self.forcefield.name)



class MSMGroup(Base):
    """
    database model describing a set of MarkovModels built during concurrently
    that share a joint clustering
    """
    __tablename__ = 'msm_groups'
    
    id = Column(Integer, primary_key=True)
    
    generators_fn = Column(String(STRING_LEN))
    microstate_selection_weights = Column(PickleType)

    def __repr__(self):
        return "<BuildRound(id={})>".format(self.id)
        