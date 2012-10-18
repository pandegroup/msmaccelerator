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
import matplotlib.pyplot as pp
import IPython as ip
import numpy as np
from msmbuilder import io
from msmbuilder import arglib
import msmbuilder.Trajectory
import shutil

from msmaccelerator import Project
from msmaccelerator.database import Session
from msmaccelerator.models import Trajectory, MSMGroup, Forcefield, MarkovModel
from msmaccelerator.utils import load_file, save_file

S = Session
pp.ion()

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

def plot_n_states(project):
    groups = S.query(MSMGroup).all()
    ids = [g.id for g in groups]
    states = [g.n_states for g in groups]
    
    pp.plot(ids, states)
    pp.title('Number of state vs. build round')
    pp.xlabel('Build round')
    pp.ylabel('Number of states')
    
    ip.embed()

def performance(project):
    trajs = S.query(Trajectory).all()
    
    for i, traj in enumerate(trajs):
        if traj.submit_time is not None:
            if traj.returned_time is not None:
                delta = traj.returned_time - traj.submit_time
                simulation_frames = traj.length - 1
                
                print 'traj %3d\tframes: %d\twall time: %s' %  (i, simulation_frames, _render_td(delta))
            else:
                print 'traj %3d\tsubmitted       %s' % (i, traj.submit_time.strftime("%d %B %r"))

def _render_td(td):
    "String representation of a timedelta as hours:minutes:seconds"
    hours, remainder = divmod(td.total_seconds(), 3600)
    minutes, seconds = divmod(remainder, 60)
    return '%02d:%02d:%02d' % (int(hours), int(minutes), int(seconds))

def main():
    parser = arglib.ArgumentParser()
    parser.add_argument('project_file', help='path to project.yaml')
    subparsers = parser.add_subparsers(dest="subparser_name")
    subparsers.add_parser('shell')
    subparsers.add_parser('performance')
    subparsers.add_parser('check_sufficient')
    subparsers.add_parser('plot_n_states')
    # subparsers.add_parser('drop_project')
    # this is only for MySQL databases which are currently not in use


    args = parser.parse_args(print_banner=False)

    project = Project(args.project_file)
    g = globals()
    if not args.subparser_name in g:
        raise NotImplementedError
    g[args.subparser_name](project)
if __name__ == '__main__':
    main()
