import os
from sqlalchemy import Column

import os
import scipy.io
import numpy as np
import msmbuilder.Trajectory
import msmbuilder.Serializer
import types


def file_store_metaclass(name, bases, dct):
    """ This metaclass injects properties into a SQLAlchemy models.
    
    It looks for sql columns named such that they start with a single
    underscore and end with "_fn", such as "_counts_fn", "_assignment_fn". For
    each of these matched columns, it also looks for a method like
    "__default_counts_fn", or "__default_assignments_fn", which should return
    a string giving the default value for that column.
    
    Then, it creates getter and setter properties, so that you can do
    >> model.counts = scipy.sparse.eye(10,10)
    or
    >> print model.counts
    <4x4 sparse matrix of type '<type 'numpy.float64'>'
	    with 4 stored elements in COOrdinate format>
    
    The formats for saving/loading that are supported are
    msmbuilder.Serializer.SaveData/msmbuilder.Serializer.LoadData,
    msmbuilder.Trajectory.LoadTrajectoryFile/msmbuilder.Serializer.LoadData,
    
    scipy.io.mmwrite/scipy.io.mmread,
    np.savetxt/np.loadtxt
    pickle.load/pickle.dump
    
    The load/save method is infered at loadtime/savetime from the extension
    of the _counts_fn or _assignments_fn database field.
    
    extesion  |  protocol
    ----------|--------
     .h5      |  msmbuilder.Serializer.LoadData/Serializer.SaveData
     .mtx     |  scipy.io.mmread/scipy.io.mmwrite
     .dat     |  np.savetxt/np.loadtxt
     .pickl   |  pickle.load/pickle.dump
     .lh5     |  msmbuilder.Trajectory.LoadTrajectoryFile
     .xtc     |                 "
     .dcd     |                 "
     .pdb     |                 "
     .txt     |  Save ascii to flat tex file
    """
    
    def save(path, value):
        ext = os.path.splitext(path)[1]

        # make directory path
        try:
            os.makedirs(os.path.split(path)[0])
        except:
            pass
            
        if os.path.exists(path):
            raise IOError(("At this point, for safety, we don't support "
                "overwriting files that already exist on disk"))
            
            
        if ext == '.h5':
            msmbuilder.Serializer.SaveData(path, value)
        elif ext == '.mtx':
            scipy.io.mmwrite(path, value)
        elif ext == '.dat':
            np.savetxt(path, value)
        elif ext == '.lh5':
            value.SaveToLHDF(path)
        elif ext == '.pdb':
            value.SaveToPDB(path)
        elif ext == '.txt':
            with open(path, 'w') as f:
                print >> f, value
        elif ext == '.xtc':
            raise NotImplementedError('I dont save xtcs')
        elif ext == '.dcd':
            raise NotImplementedError('I dont save dcds')
        else:
            raise ValueError('Extension ({}) not recognized'.format(ext))
            
    def load(path):
        ext = os.path.ext(path)[1]
        if ext == '.h5':
            return msmbuilder.Serializer.LoadData(path)
        elif ext == '.mtx':
            return scipy.io.mmread(path)
        elif ext == '.dat':
            return np.loadtxt(path)
        elif ext == '.txt':
            return open(path).read()
        elif ext in ['.lh5', '.xtc', '.dcd', 'pdb']:
            return msmbuilder.Trajectory.LoadTrajectoryFile(path)
        else:
            raise ValueError('Extension ({}) not recognized'.format(ext))
            
                
    def make_property(key, dct):
        
        # example: key=_counts_fn
        # default_filename_method = _default_counts_fn_
        
        # make sure that an attribute like _default_counts_fn_ exists
        default_filename_method = '_default' + key + '_'
        if not ((default_filename_method in dct) and \
            isinstance(dct[default_filename_method], types.FunctionType)):
            raise RuntimeError("method {} must be defined as a method".format(default_filename_method))
        
        def getter(self):
            if getattr(self, key) is None:
                return None
            print key, getattr(self, key)
            return scipy.io.mmread(getattr(self, key))
            
        def setter(self, value):
            if getattr(self, key) is None:
                default_filename = getattr(self, default_filename_method)()
                setattr(self, key, default_filename)
                
            print 'saving to ' + getattr(self, key)
            
            path = getattr(self, key)
            save(path, value)
            
            
        return property(getter, setter)
    
    
    newdct = dct.copy()   
    for k, v in dct.items():
        # look for sqlalchemy columns that start with "_" and end with "_fn"
        if k.startswith('_') and k.endswith('_fn') and isinstance(v, Column):
            basename = k[1:-3] # chop off the end characters
            newdct[basename] = make_property(k, dct)
    
    return type(name, bases, newdct)