import sys
import IPython as ip
import numpy as np
import argparse
from msmbuilder import Serializer
import msmbuilder.Trajectory

from msmaccelerator import Project
from msmaccelerator.database import Session
from msmaccelerator.models import Trajectory, MSMGroup, Forcefield, MarkovModel

S = Session

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-project_file', help='path to project.yaml')
    subparsers = parser.add_subparsers(dest="subparser_name")
    shell_parser = subparsers.add_parser('shell')
    performance_parser = subparsers.add_parser('performance')
    performance_parser = subparsers.add_parser('check_sufficient')

    args = parser.parse_args()

    project = Project(args.project_file)
    g = globals()
    if not args.subparser_name in g:
        raise NotImplementedError
    g[args.subparser_name](project)
        
    
    
def shell(project):
    print "\033[95m>> from msmaccelerator.database import Session"
    print ">> from msmaccelerator.models import (Trajectory, MSMGroup,"
    print ">>                                    Forcefield, MarkovModel)"
    print ">> P = Project({})".format(sys.argv[1])
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
    b.is_sufficient_new_data()
    
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