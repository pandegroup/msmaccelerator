# DECLARE MODELS
import os
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import types
from sqlalchemy import (Column, Integer, String, DateTime,
                        Float, ForeignKey, Table, PickleType)
Base = declarative_base()
STRING_LEN = 500

import msmbuilder.Trajectory
from Project2 import Project

class ASCII(types.TypeDecorator):
    impl = types.String
    
    def __init__(self):
        super(ASCII, self).__init__(length=500)
        
    def process_result_value(self, value, dialect):
        if isinstance(str, unicode):
            return str(value)
        return value

class Forcefield(Base):
    """
    database model describing a forcefield and how to run simulations using
    it (driver, threads, etc)
    """
    __tablename__ = 'forcefields'
    
    id = Column(Integer, primary_key=True)
    name = Column(ASCII, unique=True)
    water = Column(ASCII)
    driver = Column(ASCII)
    threads = Column(Integer)
    output_extension = Column(ASCII)
        
    def __repr__(self):
        return "<Forcefield(name={}, water={}, driver={})>".format(self.name,
            self.water, self.driver)

class Trajectory(Base):
    """
    database model describing a trajectory
    """
    
    __tablename__ = 'trajectories'
    
    id = Column(Integer, primary_key=True)
    
    init_pdb_fn = Column(ASCII)
    wqlog_fn = Column(ASCII)
    dry_xtc_fn = Column(ASCII)
    wet_xtc_fn = Column(ASCII)
    lh5_fn = Column(ASCII)
    last_wet_snapshot_fn = Column(ASCII)
    
    
    def __basefn(self):
        return '{base}/{ff}/'.format(base=Project.instance.project_dir, ff=self.forcefield.name)
    def default_init_pdb_fn(self):
       return self.__basefn() + 'Inits/{id}.pdb'.format(id=self.id)
    def default_wqlog_fn(self):
        return self.__basefn() + 'WQLogs/{id}.log.txt'.format(id=self.id)
    def default_dry_xtc_fn(self):
        return self.__basefn() + 'DryXTCs/{id}.xtc'.format(id=self.id)
    def default_wet_xtc_fn(self):
        return self.__basefn() + 'WetXTCs/{id}.xtc'.format(id=self.id)
    def default_lh5_fn(self):
        return self.__basefn() + 'LH5s/{id}.lh5'.format(id=self.id)
    def default_last_wet_snapshot_fn(self):
        return self.__basefn() + 'LastWetSnapshots/{id}.lh5'.format(id=self.id)

    def populate_default_filenames(self):
        if self.id is None:
            raise Exception(('self.id is None!. Did you forget to commit before '
                'calling this method?'))
                
        for name in ['init_pdb_fn', 'wqlog_fn', 'dry_xtc_fn', 'wet_xtc_fn',
            'lh5_fn', 'last_wet_snapshot_fn']:
            if getattr(self, name) is None:
                # call the method
                default = getattr(self, 'default_' + name)()
                # set the field
                setattr(self, name, default)
                
                #just to be safe, lets make the directories if they dont exist
                try:
                    os.makedirs(os.path.split(default)[0])
                except OSError:
                    pass
                
                
        
    host = Column(ASCII)
    mode = Column(ASCII)
    name = Column(ASCII)
    returned_time = Column(Float)
    submit_time = Column(Float)
    length = Column(Integer) # in frames
    
    # many to one relationship
    forcefield_id = Column(Integer, ForeignKey('forcefields.id'))
    forcefield = relationship("Forcefield", backref=backref('trajectories', order_by=id))
    
    def __repr__(self):
        return "<Trajectory(name={}, init_pdb_fn={})>".format(self.name,
            self.init_pdb_fn)
        
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
            conf = msmbuilder.Trajectory.LoadTrajectoryFile(self.last_wet_snapshot_fn)
            xyz = msmbuilder.Trajectory.ReadFrame(self.wet_xtc_fn, frame_index)
        else:
            conf = msmbuilder.Trajectory.LoadTrajectoryFile(Project().pdb_topology_file)
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
    
    id = Column(Integer, primary_key=True)
    counts_fn = Column(ASCII)
    assignments_fn = Column(ASCII)
    inverse_assignments_fn = Column(ASCII)
    
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
    
    generators_fn = Column(ASCII)
    microstate_selection_weights = Column(PickleType)

    def __repr__(self):
        return "<BuildRound(id={})>".format(self.id)
        