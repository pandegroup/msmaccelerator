import yaml
import numpy as np
from pprint import pprint
from Emsmbuilder import metrics
from Emsmbuilder import clustering
from functools import partial
import IPython as ip

def parse_yaml_input(input_file):
    # Set up the yaml parser
    def array_constructor(loader, node):
        value = loader.construct_scalar(node)
        path, dtype = map(unicode.strip, value.split(','))
        return np.loadtxt(path, dtype)
    def metric_loader(loader, node):
        value = loader.construct_scalar(node)
        attr = getattr(metrics, value)
        return attr
    def clustering_loader(loader, node):
        value = loader.construct_scalar(node)
        attr = getattr(clustering, value)
        return attr
    class reprd_partial(partial):
        def __repr__(self):
            args = ', '.join(repr(a) for a in self.args) if len(self.args) != 0 else None
            kwargs = ', '.join('%s=%s' % (repr(k), repr(v))  for k, v in self.keywords.iteritems())
            together = ', '.join(filter(lambda i: i is not None, (args, kwargs)))
            return "functools.partial(%r, %s)" % (self.func, together)
    
    yaml.add_constructor('!array', array_constructor)
    yaml.add_constructor('!metric', metric_loader)
    yaml.add_constructor('!clustering', clustering_loader)
    
    # parse the yaml file
    params = yaml.load(input_file)
    
    # construct the metric and clustering
    metric = params['metric']['type'](**params['metric']['init_kwargs'])
    params['clustering'] = reprd_partial(params['clustering']['type'], metric=metric, **params['clustering']['init_kwargs'])
    # params['clustering] is now a slightly confusing entity. It's a partially
    # applied constructor to the clustering object's constructor. To use it
    # you simply provide the trajectories like this (you need to use the keyword)
    # params['clustering'](trajectories=my_list_of_trajectories) which returns
    # a clusterer object (Emsmbuilder.clustering.KCenters, etc) which you can
    # call get_generators_as_traj() on or something
        
    return params
    
params = parse_yaml_input(open('input_file.yaml'))
pprint(params)
