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
from threading import Lock
Session = scoped_session(sessionmaker())


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


_session_lock = Lock()
def with_db_lock(f):
    @functools.wraps(f)
    def wrap(*args, **kwargs):
        _session_lock.acquire()
        try:
            return f(*args, **kwargs)
        finally:
            _session_lock.release()
    return wrap