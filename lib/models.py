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

# DECLARE MODELS
import os
import numpy as np
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlalchemy import types
from sqlalchemy import (Column, Integer, String, DateTime,
                        Float, ForeignKey, Table, PickleType,
                        Boolean)
Base = declarative_base()
STRING_LEN = 500

import msmaccelerator.Project
import msmbuilder.Trajectory

class _PopulateMixin(object):
    """Mixin superclass that adds the populate_default_filenames method.
    
    The idea here is to make it easier to store files in the database by only
    storing the path to the filename on disk.
    
    models that inherit from _PopulateMixin should have string columns such as
    lh5_fn or counts_fn (ending in "_fn"), and a corresponding method called
    default_lh5_fn or default_counts_fn. When populate_default_filenames is
    called, those filenames will be set by calling the default function.
    
    This gives you a little bit more flexibility than what you can accomplish using
    SQLAlchemey defaults alone, because you choose when populate_default_filenames
    is called.
    
    This enables you to set the defaults only *AFTER* the model has been flushed
    into the database, assigning it an ID. This makes it possible for the default
    value of one of the filenames to depend on the id or on some other relationship.

    I know this is not a very elegant solution, but I'm not sure how to get a
    truely transparent hybrid of in-database and on-disk storage with all of the
    aspects like default filenames that depend on id or relationships...
    """
    
    def populate_default_filenames(self):
        if self.id is None:
            raise Exception(('self.id is None!. Did you forget to commit before '
                'calling this method?'))
        
        # look for Columns that end with "_fn"; there should be
        # a corresponding default_whatever_fn method
        for name, val in self.__class__.__dict__.items():
            if name.endswith('_fn') and isinstance(val, InstrumentedAttribute):
                # call the method
                if not hasattr(self, 'default_' + name):
                    raise ValueError()
                default = getattr(self, 'default_' + name)()
                # set the field
                setattr(self, name, default)
                
                #just to be safe, lets make the directories if they dont exist
                try:
                    os.makedirs(os.path.split(default)[0])
                except OSError:
                    pass
    

class ASCII(types.TypeDecorator):
    """
    Honestly, I don't know if this actually works. Its supposed to be a
    string database type that ensures that the values coming out of the db
    are of type str (not of type unicode), but it doesn't seem to work.
    
    I haven't checked closely whats going on (RTM 8/6)
    """
    impl = types.String
    
    def __init__(self):
        super(ASCII, self).__init__(length=500)
        
    def process_result_value(self, value, dialect):
        if isinstance(str, unicode):
            return str(value)
        return value

class Forcefield(Base, _PopulateMixin):
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
    cost = Column(Float)
    output_extension = Column(ASCII)
    true_kinetics = Column(Boolean)
        
    def __repr__(self):
        return "<Forcefield(name={}, water={}, driver={})>".format(self.name,
            self.water, self.driver)

class Trajectory(Base, _PopulateMixin):
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
        return '{base}/{ff}/'.format(base=msmaccelerator.Project.instance.project_dir, ff=self.forcefield.name)
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
        return self.__basefn() + 'LastWetSnapshots/{id}.pdb'.format(id=self.id)

    host = Column(ASCII)
    mode = Column(ASCII)
    name = Column(ASCII)
    returned_time = Column(DateTime)
    submit_time = Column(DateTime)
    length = Column(Integer) # in frames

    initialized_from_id = Column(Integer)
    initialized_from_frame = Column(Integer)
    
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
            conf = msmbuilder.Trajectory.load_trajectory_file(str(self.last_wet_snapshot_fn))
            xyz = msmbuilder.Trajectory.read_frame(str(self.wet_xtc_fn), frame_index)
        else:
            conf = msmbuilder.Trajectory.load_trajectory_file(Project().pdb_topology_file)
            xyz = msmbuilder.Trajectory.read_frame(str(xtc_path), frame_index)

        if not (xyz.shape == conf['XYZList'][0,:,:].shape):
            raise Exception("Number of atoms is wrong: xyz.shape={0}, conf['XYZList'][0,:,:].shape={1}".format(xyz.shape, conf['XYZList'][0,:,:].shape))

        conf['XYZList'] = np.array([xyz])

        return conf



class MarkovModel(Base, _PopulateMixin):
    """
    database model describing an MSM built from a set of trajectories, all in the
    same forcefield
    
    Each BuildRound contains multiple MarkovModels
    """
    __tablename__ = 'markov_models'
    
    id = Column(Integer, primary_key=True)
    counts_fn = Column(ASCII)
    assignments_fn = Column(ASCII)
    distances_fn = Column(ASCII)
    inverse_assignments_fn = Column(ASCII)
    
    # I'm putting the id in the default filename to guard against of off chance
    # that two MarkovModels with the same forcefied in the same round get created
    # and then they would otherwise be given the same filenames
    def __basefn(self):
        return '{base}/models/{ff}/{round}/'.format(base=msmaccelerator.Project.instance.project_dir,
            ff=self.forcefield.name, round=self.msm_group_id)
    def default_counts_fn(self):
        return self.__basefn() + 'tCounts.{id}.mtx'.format(id=self.id)
    def default_assignments_fn(self):
        return self.__basefn() + 'Assignments.{id}.h5'.format(id=self.id)
    def default_distances_fn(self):
        return self.__basefn() + 'Distances.{id}.h5'.format(id=self.id)
    def default_inverse_assignments_fn(self):
        return self.__basefn() + 'InvAssignmnets.{id}.pickl'.format(id=self.id)
    
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
    microstate_selection_weights = Column(PickleType)
    
    def __repr__(self):
        return "<MarkovModel(id={}, msm_group={}, forcefield={})>".format(self.id,
            self.msm_group_id, self.forcefield.name)



class MSMGroup(Base, _PopulateMixin):
    """
    database model describing a set of MarkovModels built during concurrently
    that share a joint clustering
    """
    __tablename__ = 'msm_groups'
    
    id = Column(Integer, primary_key=True)
    generators_fn = Column(ASCII)
    n_states = Column(Integer)
    
    def __basefn(self):
        return '{base}/models/gens/'.format(base=msmaccelerator.Project.instance.project_dir)
    def default_generators_fn(self):
        return self.__basefn() + 'Gens.{id}.lh5'.format(id=self.id)
    
    trajectories = relationship("Trajectory",  backref='msm_groups',
        secondary=Table('msm_group_trajectories', Base.metadata,
                Column('markov_model_id', Integer, ForeignKey('msm_groups.id')),
                Column('trajectory_id', Integer, ForeignKey('trajectories.id'))))
    

    def __repr__(self):
        return "<BuildRound(id={})>".format(self.id)