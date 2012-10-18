#!/usr/bin/env python
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



import sys
import IPython as ip
import numpy as np
#from msmbuilder import io
from msmbuilder import arglib
import msmbuilder.Trajectory
import shutil

from msmaccelerator import Project
from msmaccelerator.database import Session
from msmaccelerator.models import Trajectory, MSMGroup, Forcefield, MarkovModel
from msmaccelerator.utils import load_file, save_file

S = Session

def main():
    parser = arglib.ArgumentParser()
    parser.add_argument('project_file', help='path to project.yaml')
    subparsers = parser.add_subparsers(dest="subparser_name")
    subparsers.add_parser('shell')
    subparsers.add_parser('performance')
    subparsers.add_parser('check_sufficient')
    subparsers.add_parser('drop_project')
    subparsers.add_parser('cleanup')


    args = parser.parse_args(print_banner=False)

    project = Project(args.project_file)
    g = globals()
    if not args.subparser_name in g:
        raise NotImplementedError
    g[args.subparser_name](project)
        
    
def drop_project(project):
    val = raw_input('Are you sure you want to delete the project (deleting data and db)? yes/[no] ')
    if val != 'yes':
        print 'Exiting.'
        return

    connection = S.connection()
    connection.execute('DROP DATABASE {};'.format(project.mysql_db))
    print 'Dropped database.'
    shutil.rmtree(project.project_dir)
    print 'Files deleted'
    

def shell(project):
    print "\033[95m>> from msmaccelerator.database import Session"
    print ">> from msmaccelerator.models import (Trajectory, MSMGroup,"
    print ">>                                    Forcefield, MarkovModel)"
    print ">> P = Project(<your_project_file>)"
    print ">> S = Session"
    print "\033[0m"
    print "\n"*2
    
    P = project
    from IPython.frontend.terminal.embed import InteractiveShellEmbed
    ipshell = InteractiveShellEmbed(banner1='')
    
    ipshell()
    
def check_sufficient(project):
    from msmaccelerator import Builder
    b = Builder(project)
    print 'Running\n'
    b.is_sufficient_new_data()
    print 'Done'
    
def performance(project_file):
    trajs = S.query(Trajectory).all()
    
    
    for i, traj in enumerate(trajs):
        if traj.submit_time is not None:
            if traj.returned_time is not None:
                delta = traj.returned_time - traj.submit_time
                simulation_frames = traj.length - 1
                
                print 'traj %3d\tframes: %d\twall time:%s' %  (i, simulation_frames, str(delta))
            else:
                print 'traj %3d\tsubmitted       %s' % (i, traj.submit_time.strftime("%d %B %r"))

    

if __name__ == '__main__':
    main()
