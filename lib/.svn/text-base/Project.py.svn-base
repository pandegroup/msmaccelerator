"""
Encapsulate some of the file handling stuff
"""
import yaml
import os, sys
import functools
import numpy as np
from glob import glob
import re
import threading
from msmbuilder import metrics
from msmbuilder import clustering
from msmbuilder.Trajectory import Trajectory
from monakos.clustering import MonakosHierarchical
from monakos import Sampling
clustering.MonakosHierarchical = MonakosHierarchical

class Project(object):
    def __init__(self, params_fn):
        """Load up by parse the yaml parameter file"""

        self.params_fn = params_fn # save this incase we need it --TJL
        
        # Set up the yaml parser
        def array_constructor(loader, node):
            value = loader.construct_scalar(node)
            path, dtype = map(unicode.strip, value.split(','))
            path = os.path.join(params_dir, path)
            return np.loadtxt(path, dtype)
        def metric_loader(loader, node):
            value = loader.construct_scalar(node)
            attr = getattr(metrics, value)
            return attr
        def clustering_loader(loader, node):
            value = loader.construct_scalar(node)
            attr = getattr(clustering, value)
            return attr
        class reprd_partial(functools.partial):
            def __repr__(self):
                args = ', '.join(repr(a) for a in self.args) if len(self.args) != 0 else None
                kwargs = ', '.join('%s=%s' % (repr(k), repr(v))  for k, v in self.keywords.iteritems())
                together = ', '.join(filter(lambda i: i is not None, (args, kwargs)))
                return "functools.partial(%r, %s)" % (self.func, together)
        
        yaml.add_constructor('!array', array_constructor)
        yaml.add_constructor('!metric', metric_loader)
        yaml.add_constructor('!clusterer', clustering_loader)
        
        # parse the yaml file
        params_dir = os.path.abspath(os.path.split(params_fn)[0])
        params = open_yaml(params_fn)
        
        self.metric = params['metric']['type'](**params['metric']['init_kwargs'])
        self.clusterer = reprd_partial(params['clusterer']['type'], metric=self.metric, **params['clusterer']['init_kwargs'])
        
        # self.clustering is now a slightly confusing entity. It's a partially
        # applied constructor to the clustering object's constructor. To use it
        # you simply provide the trajectories like this (you need to use the keyword)
        # project.clusterer(trajectories=my_list_of_trajectories) which returns
        # a clusterer object (msmbuilder.clustering.KCenters, etc) which you can
        # call get_generators_as_traj() on or something
        
        try:
            self.forcefields = params['forcefields']
            for ff in self.forcefields:
                ff['driver'] = os.path.join(params_dir, os.path.expanduser(ff['driver']))
            
            self.pdb_topology_file = os.path.join(params_dir, params['pdb_topology_file'])
            
            # if there is a parameter called starting_confs_lh5, we'll randomly pull
            # conformations from there to start the initial equilibration runs from
            # otherwise, we'll just use the pdb_topology_file
            if 'starting_confs_lh5' in params:
                self.starting_confs_lh5 = os.path.join(params_dir, params['starting_confs_lh5'])
            else:
                self.starting_confs_lh5 = None
            
            self.method = params['method']
            assert self.method in ['implicit', 'explicit']
            
            self.stride = params['stride']
            self.symmetrize = params['symmetrize']
            self.lagtime = params['lagtime']
            self.trim = params['trim']
            self.cross_assign = params['cross_assign']
            self.adaptive_sampling = params['adaptive_sampling']
            self.num_rounds = params['num_rounds']
            self.rootdir = params['project_dir']
            self.num_trajs_sufficient_for_round = params['num_trajs_sufficient_for_round']
        except KeyError as e:
            raise KeyError('A required option, %s, was missing from params file' % e)
        
        if 'xtc_freq_ps' in params:
            raise Exception('xtc_freq_ps is deprecated.')
        
        self.rootdir = os.path.join(params_dir, self.rootdir)
        if not os.path.exists(self.rootdir):
            os.mkdir(self.rootdir)
            
        self.builder_log_fn = os.path.join(self.rootdir, 'BuilderLog.yaml')
        self.traj_log_fn = os.path.join(self.rootdir, 'TrajLog.yaml')
        if not os.path.exists(self.builder_log_fn):
            write_yaml([], self.builder_log_fn)
        if not os.path.exists(self.traj_log_fn):
            write_yaml([], self.traj_log_fn)
        
        self._builder_log_lock = threading.Lock()
        self._traj_log_lock = threading.Lock()
        
        if len(self.forcefields) > 1:
            self.joint_ff = {'name': 'JointFF'}
            
        self.validate_params()
        
        
    
    def validate_params(self):
        Trajectory.LoadTrajectoryFile(self.pdb_topology_file)
        
        e0 = 'stride must be an integer greater than or equal to one'
        assert self.stride == int(self.stride), e0
        assert 1 <= self.stride, e0
            
        if self.symmetrize is not None:
            self.symmetrize = str(self.symmetrize).lower() 
        else:
            self.symmetrize = 'none'
        assert self.symmetrize in ['mle', 'transpose', 'none'], ('Symmetrize must '
            'be one of "mle", "transpose", or "none". Refer to the MSMbuilder ' 
            'docs for details')
        
        e1 = 'lagtime must be an integer greater than or equal to 1'
        assert self.lagtime == int(self.lagtime), e1
        assert 1 <= self.lagtime, e1
        
        e2 = 'cross_assign must be a boolean'
        assert self.cross_assign == bool(self.cross_assign), e2
        
        e3 = 'adaptive_sampling must be the name (string) of a method in'
        assert self.adaptive_sampling == str(self.adaptive_sampling)
        try:
            getattr(Sampling, self.adaptive_sampling)
        except:
            raise AssertionError(e3)
                
        e5 = 'num rounds must be an integer greater than or equal to zero'
        assert self.num_rounds == int(self.num_rounds), e5
        assert self.num_rounds >= 0, e5
        
        num_gold_standards = 0
        for ff in self.forcefields:
            assert ff['name'] != 'AllFF', 'no forcefield can be named "AllFF"'
            assert os.path.exists(ff['driver']), "can't find driver"
            assert 'output_extension' in ff, 'need to supply output format'
            assert ff['output_extension'] in ['.xtc', '.dcd'], 'output format needs to be .xtc or .dcd'
            assert 'threads' in ff, 'ff must have a field called threads'
            assert 'water' in ff, 'ff must have field called water'
            
            if 'gold_standard' in ff and ff['gold_standard'] is True:
                num_gold_standards += 1
            else:
                ff['gold_standard'] = False
            
            e7 = 'each ff needs a cost, which needs to be an int'
            assert 'cost' in ff and ff['cost'] == int(ff['cost']), e7
            
        
        assert num_gold_standards == 1, '1 and only 1 forcefield needs to be labeled gold_standard'
        
        assert self.num_trajs_sufficient_for_round == int(self.num_trajs_sufficient_for_round), 'sdfsdsdsdfsdf'
            
    
    def ff_name_2_index(self, ff_name):
        for i, ff in enumerate(self.forcefields):
            if ff['name'] == ff_name:
                return i
        raise IndexError('Not Found')

    @property
    def builder_log(self):
        with self._builder_log_lock:
            val = open_yaml(self.builder_log_fn)
        return val
        
    @property
    def traj_log(self):
        with self._traj_log_lock:
            val = open_yaml(self.traj_log_fn)
        return val
    
    def builder_log_add(self, round_dict):
        log = self.builder_log
        with self._builder_log_lock:
            log.append(round_dict)
            write_yaml(log, self.builder_log_fn)

    def traj_log_add(self, traj_dict):
        log = self.traj_log
        with self._traj_log_lock:
            log.append(traj_dict)
            write_yaml(log, self.traj_log_fn)

    @property
    def rounds_on_disk(self):
        rounds = [[] for ff in self.forcefields]
        for i, ff in enumerate(self.forcefields):
            ps = glob(os.path.join(self.rootdir, ff['name'], 'round*'))
            if len(rounds) > 0:
                rounds[i] = [int(re.search('round(\d+)', p).group(1)) for p in ps]
                rounds[i].sort()
                
        return rounds
    
    def starting_conf_dir(self, ff_name):
        return self._ff_dir(ff_name, 'Inits')
    
    def last_snapshot_dir(self, ff_name):
        return self._ff_dir(ff_name, 'LastSnapshots')

    def xtc_dir(self, ff_name):
        return self._ff_dir(ff_name, 'XTCs')

    def wqlog_dir(self, ff_name):
        return self._ff_dir(ff_name, 'WQLogs')

    def traj_dir(self, ff_name):
        return self._ff_dir(ff_name, 'Trajectories')

    def _ff_dir(self, ff_name, name):
        val = os.path.join(self.rootdir, ff_name, name)
        if not os.path.exists(val):
            os.makedirs(val)
        return val

    def _ff_round_fn(self, ff_name, round_num, fn):
        val = os.path.join(self.rootdir, ff_name, 'round%d' % round_num, fn)
        directory = os.path.split(val)[0]
        if not os.path.exists(directory):
            #print 'making directory: %s' % directory
            os.makedirs(directory)
        return val

    def counts_fn(self, ff_name, round_num):
        return self._ff_round_fn(ff_name, round_num, 'tCounts.mtx')

    def assignments_fn(self, ff_name, round_num):
        return self._ff_round_fn(ff_name, round_num, 'Assignments.h5')

    def medoids_fn(self, ff_name, round_num):
        return self._ff_round_fn(ff_name, round_num, 'Gens.lh5')

def open_yaml(fn):
    #try:
    with open(fn) as f:
        a = yaml.load(f)
    #except Exception as e:
    #    raise IOError('Something went wrong with yaml: {0}'.format(e.message))
    return a
    
def write_yaml(obj, fn):
    with open(fn, 'w') as f:
        yaml.dump(obj, f)
