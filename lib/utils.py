import functools
import collections
from itertools import ifilterfalse
from operator import itemgetter
import logging
import logging.handlers
import msmbuilder.Trajectory
import msmbuilder.Serializer
import scipy.io
import cPickle as pickle
import sys, os

def load_file(path):
    """Load a file
    
    .mtx, .h5, .lh5, .pdb, .pickle
    """
    path = str(path)
    
    ext = os.path.splitext(path)[1].lower()
    if ext == '.mtx':
        return scipy.io.mmread(path)
    elif ext == '.h5':
        return msmbuilder.Serializer.LoadData(path)
    elif ext in ['.lh5', '.pdb']:
        return msmbuilder.Trajectory.LoadTrajectoryFile(path)
    elif ext == '.pickl':
        return pickle.load(open(path))
    else:
        raise NotImplementedError(ext)

def save_file(path, value):
    """Save a file
    
    .mtx, .h5, .lh5, .xtc, .pdb, .pickl
    """
    path = str(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == '.mtx':
        return scipy.io.mmwrite(path, value)
    elif ext == '.h5':
        return msmbuilder.Serializer.SaveData(path, value)
    elif ext == '.lh5':
        return value.SaveToLHDF(path)
    elif ext == '.xtc':
        return value.SaveToXTC(path)
    elif ext == '.pdb':
        return value.SaveToPDB(path)
    elif ext == '.pickl':
        return pickle.dump(value, open(path, 'w'))
    else:
        raise NotImplementedError(ext)


#http://stackoverflow.com/questions/1363839/python-singleton-object-instantiation
class Singleton(type):
    def __init__(self, name, bases, dict):
        super(Singleton, self).__init__(name, bases, dict)
        self.instance = None

    def __call__(self, *args, **kw):
        if self.instance is None:
            self.instance = super(Singleton, self).__call__(*args, **kw)

        return self.instance
        
        


class GMailHandler(logging.handlers.SMTPHandler):
    "Logging handler to send email from an msmaccelerator account"
    def __init__(self, toaddrs, subject):
        mailhost = 'smtp.gmail.com'
        mailport = 587
        gmail_user = 'msmaccelerator@gmail.com'
        gmail_password = 'thisisthepassword'
        empty_tuple = ()

        logging.handlers.SMTPHandler.__init__(self, (mailhost, mailport),
            gmail_user, toaddrs, subject, (gmail_user, gmail_password), empty_tuple)

class Counter(dict):
    'Mapping where default values are zero'
    def __missing__(self, key):
        return 0

def lru_cache(maxsize=100):
    '''Least-recently-used cache decorator.

    Arguments to the cached function must be hashable.
    Cache performance statistics stored in f.hits and f.misses.
    Clear the cache with f.clear().
    http://en.wikipedia.org/wiki/Cache_algorithms#Least_Recently_Used

    '''
    maxqueue = maxsize * 10
    def decorating_function(user_function,
            len=len, iter=iter, tuple=tuple, sorted=sorted, KeyError=KeyError):
        cache = {}                  # mapping of args to results
        queue = collections.deque() # order that keys have been used
        refcount = Counter()        # times each key is in the queue
        sentinel = object()         # marker for looping around the queue
        kwd_mark = object()         # separate positional and keyword args

        # lookup optimizations (ugly but fast)
        queue_append, queue_popleft = queue.append, queue.popleft
        queue_appendleft, queue_pop = queue.appendleft, queue.pop

        @functools.wraps(user_function)
        def wrapper(*args, **kwds):
            # cache key records both positional and keyword args
            key = args
            if kwds:
                key += (kwd_mark,) + tuple(sorted(kwds.items()))

            # record recent use of this key
            queue_append(key)
            refcount[key] += 1

            # get cache entry or compute if not found
            try:
                result = cache[key]
                wrapper.hits += 1
            except KeyError:
                result = user_function(*args, **kwds)
                cache[key] = result
                wrapper.misses += 1

                # purge least recently used cache entry
                if len(cache) > maxsize:
                    key = queue_popleft()
                    refcount[key] -= 1
                    while refcount[key]:
                        key = queue_popleft()
                        refcount[key] -= 1
                    del cache[key], refcount[key]

            # periodically compact the queue by eliminating duplicate keys
            # while preserving order of most recent access
            if len(queue) > maxqueue:
                refcount.clear()
                queue_appendleft(sentinel)
                for key in ifilterfalse(refcount.__contains__,
                                        iter(queue_pop, sentinel)):
                    queue_appendleft(key)
                    refcount[key] = 1


            return result

        def clear():
            cache.clear()
            queue.clear()
            refcount.clear()
            wrapper.hits = wrapper.misses = 0

        wrapper.hits = wrapper.misses = 0
        wrapper.clear = clear
        return wrapper
    return decorating_function


def setup_logger(file_level, email_level, filename, email_addrs, email_subject='MSMAccelerator Log'):
    logger = logging.getLogger('MSMAccelerator')
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(name)s: %(asctime)s: %(levelname)s: %(message)s', '%m/%d/%y [%H:%M:%S]')

    if len(email_addrs) > 0:
        gmail_handler = GMailHandler(email_addrs, email_subject)
        gmail_handler.setLevel(email_level)
        gmail_handler.setFormatter(formatter)
        logger.addHandler(gmail_handler)

    file_handler = logging.FileHandler(filename)
    file_handler.setLevel(file_level)
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(file_level)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)

    return logger