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

import models
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import create_engine
import sqlalchemy.exc
import functools
import threading
Session = scoped_session(sessionmaker(autoflush=False))


def _connect_to_mysql_db(user, password, host, db):
    """
    Connect to a mysql database
    """
    
    db_path = "mysql://{}:{}@{}/{}".format(user, password, host, db)
    db_path2 = "mysql://{}:{}@{}".format(user, password, host)
    
    def connect():
        engine = create_engine(db_path, echo=False)
        Session.configure(bind=engine)
        models.Base.metadata.create_all(engine)
    def create():
        engine = create_engine(db_path2, echo=True)
        connection = engine.connect()
        connection.execute('CREATE DATABASE %s' % db)
        connection.close()
    
    try:
        connect()
    except sqlalchemy.exc.OperationalError:
        create()
        connect()


def _connect_to_sqlite_db(db_path):
    #db_path =  os.path.join(self.project_dir, 'db.sqlite')
    engine = create_engine('sqlite:///{}'.format(db_path), echo=False)
    Session.configure(bind=engine)
    models.Base.metadata.create_all(engine)


# @with_db_lock decorator 
# ======================== #
# this decorator is a little complicated, and is used to get around the fact
# that sqlite3 does not support concurrent access from multiple threads

# this decorator should be used by all functions that directly use the session
# object.

# first, in manages a reentrant lock associated with the database. only
# one thread can hold the lock at a time. This is actually probably not doing
# much, since we're already using a sqlalchemy "scoped session", which is supposed
# to do this already

# the second thing is that this decoractor runs a commit on the Session after
# the function exits. SORRY. this is definitly not the most IO efficient, but
# it appears to be necessary (empirically, i get a DB locking error if I don't
# have it). Basically, with SQLite, a transaction might still be open in the session
# even after the function finishes and the lock is released -- it has no idea that
# the function returned. But commit()ing seems to end the transaction.


_session_lock = threading.RLock()
def with_db_lock(f):
    @functools.wraps(f)
    def wrap(*args, **kwargs):
        #print 'Acquiring lock! %s' % threading.current_thread()
        _session_lock.acquire()
        try:
            return_value = f(*args, **kwargs)
        finally:
            #print 'Releasing Lock! %s' % threading.current_thread()
            Session.commit()
            _session_lock.release()
        return return_value
    return wrap