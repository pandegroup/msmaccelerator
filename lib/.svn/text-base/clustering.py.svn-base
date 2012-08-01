from msmbuilder.clustering import Hierarchical
class MonakosHierarchical(Hierarchical):
    def __init__(self, metric, trajectories, k=None, cutoff_distance=None, **kwargs):
        self.k = k
        self.cutoff_distance = cutoff_distance
        super(MonakosHierarchical, self).__init__(metric, trajectories, **kwargs)
    
    def get_assignments(self):
        return super(MonakosHierarchical, self).get_assignments(self.k, self.cutoff_distance)
        
        
    def get_distances(self):
        raise NotImplementedError(("Hierarhical clustering doesnt have distances to "
                "generators, since there are no generators"))