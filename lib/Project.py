"""
Encapsulate some of the file handling stuff
"""
import yaml
import os, sys
import functools
import numpy as np
from msmbuilder import metrics, clustering
from msmbuilder.Trajectory import Trajectory
from sqlalchemy import create_engine

import models
import sampling
from utils import Singleton
from database import Session


class Project(object):
    __metaclass__ = Singleton
    # SINGLETON OBJECT
    #
    # The first time that this object gets created, the __init__
    # method will run, but subsequently, calls to Project() will
    # just get you the same copy of the object back. We only want
    # one copy of the Project around at once, because we want everything
    # to have the same database connection
    
    
    def __init__(self, params_fn):
        """Load up by parse the yaml parameter file"""
        
        self.params_fn = params_fn
        self.params_dir = os.path.abspath(os.path.split(params_fn)[0])
        
        # Set up the yaml parser
        def array_constructor(loader, node):
            value = loader.construct_scalar(node)
            path, dtype = map(unicode.strip, value.split(','))
            path = os.path.join(self.params_dir, path)
            return np.loadtxt(path, dtype)
        def metric_loader(loader, node):
            value = loader.construct_scalar(node)
            attr = getattr(metrics, value)
            return attr
        def clustering_loader(loader, node):
            value = loader.construct_scalar(node)
            attr = getattr(clustering, value)
            return attr

        yaml.add_constructor('!array', array_constructor)
        yaml.add_constructor('!metric', metric_loader)
        yaml.add_constructor('!clusterer', clustering_loader)
        
        # parse the yaml file


        with open(params_fn) as f:
            params = yaml.load(f)

        
        self.metric = params['metric']['type'](**params['metric']['init_kwargs'])
        self.clusterer = functools.partial(params['clusterer']['type'], metric=self.metric, **params['clusterer']['init_kwargs'])
        
        # self.clustering is now a slightly confusing entity. It's a partially
        # applied constructor to the clustering object's constructor. To use it
        # you simply provide the trajectories like this (you need to use the keyword)
        # project.clusterer(trajectories=my_list_of_trajectories) which returns
        # a clusterer object (msmbuilder.clustering.KCenters, etc) which you can
        # call get_generators_as_traj() on or something
    
        self.pdb_topology_file = os.path.join(self.params_dir, params['pdb_topology_file'])
        if 'starting_confs_lh5' in params:
            self.starting_confs_lh5 = os.path.join(self.params_dir, params['starting_confs_lh5'])
        else:
            self.starting_confs_lh5 = None

        self.method = params['method']
        self.stride = params['stride']
        self.symmetrize = params['symmetrize']
        self.lagtime = params['lagtime']
        self.trim = params['trim']
        self.n_rounds = params['num_rounds']
        self.project_dir = params['project_dir']
        self.num_trajs_sufficient_for_round = params['num_trajs_sufficient_for_round']
        
        self.adaptive_sampling = getattr(sampling, params['adaptive_sampling'])
        
        self.__validate()
        self.__connect_to_db()
        self.add_forcefields_to_db(params['forcefields'])
        
    
    def __validate(self):
        if not self.method in ['implicit', 'explicit']:
            raise ValueError(self.method)
            
        if not os.path.exists(self.project_dir):
            os.makedirs(self.project_dir)
            
        if not self.symmetrize in ['mle', 'transpose', False, None]:
            raise ValueError(self.symmetrize)
            
            
    def __connect_to_db(self):       
        db_path =  os.path.join(self.project_dir, 'db.sqlite')
        engine = create_engine('sqlite:///{}'.format(db_path), echo=False)
        Session.configure(bind=engine)
        models.Base.metadata.create_all(engine) 
        
        
    def add_forcefields_to_db(self, p):
        if Session.query(models.Forcefield).count() == 0:
            # add forcefields
            for ff in p:
                # obj = models.Forcefield(name=ff['name'], water=ff['water'],
                #     driver=os.path.join(self.params_dir, ff['driver']),
                #     output_extension=ff['output_extension'],
                #     threads=ff['threads'], cost=ff['cost'],
                #     true_kinetics=ff['true_kinetics'])
                obj = models.Forcefield(**ff)
                    
                Session.add(obj)

            Session.commit()

        else:
            print "NOTE: I'M NOT PARSING NEW FORCEFIELDS"
        
        
# if there is a parameter called starting_confs_lh5, we'll randomly pull
# conformations from there to start the initial equilibration runs from
# otherwise, we'll just use the pdb_topology_file

# self.starting_confs_lh5 = os.path.join(params_dir, params['starting_confs_lh5'])
# self.pdb_topology_file = os.path.join(params_dir, params['pdb_topology_file'])
# self.starting_confs_lh5 = None
        
    