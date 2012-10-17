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
import os
import subprocess

"""
Prototype driver class, for AMBER. This script is the part of MSMAccelerator that 
gets shipped off to workers so they can actually execute MD runs.

This driver has been hardcoded to use the following:

-- A heating run as part of the equilibration, taking the system from 0 to 
   the simulation temperature
-- Langevin integration with gamma = 2.0 ps^{-1}
-- NPT equilibration with the Berendsen thermostat
-- NVT production

This script assumes you have an executable called "pmemd.MPI" on the cluster
used to generate MD runs. It will look for this in your $PATH, along with mpirun,
and employ those guys to get the job done.

This should allow for the following:

 $ python amber_driver.py <pdb_fn> <forcefield> <water> <mode> <threads> [xtcout_fn]

This script must be entirely self-contained to function properly. Therefore, it
contains some parameters that must be set manually. Manual review of this is a
necessary evil.

We've found that for compatability with dumb cluster sysadmins, this script should
be able to run with old-ass python (e.g., python2.4). We recommend maintaining this
compatability.
"""

# + ------------    Manually set parameters   ------------ +
#
# this constitutes a minimalist set of parameters you
# may want to use in an MD run - further control
# can be gained by directly modifying the .mdp text in the
# write_mdp method below

PRODUCTION_SETTINGS = {
    'timestep'       : 0.002,         # ps
    'total_time'     : 1,             # ps / 1ps
    'log_freq'       : 1,             # timesteps / 20ps
    'xtc_freq'       : 1,             # timesteps
    'temperature'    : 300,           # Kelvin
    'box_size'       : 3.5,           # Nanometers
    'Cl'             : 0,             # number Cl's
    'Na'             : 0              # number Na's
}

# the equilibration run copys over any parameters not
# specified from the Production set. So just add any-
# thing that should be changed for production

EQUILIBRATION_SETTINGS = {
    'timestep'       : 0.002,         # ps
    'total_time'     : 1,             # ps
    'barostat'       : 'berendsen',   # type
    'pressure'       : 1              # bar
}

# + ----------    END Manually set parameters   ---------- +


class Parameters:

    def __init__(self, forcefield, water, mode):
        self.forcefield = forcefield
        self.water      = water
        self.mode       = mode.capitalize()     # production or equilibration
        
        self.Production = PRODUCTION_SETTINGS
        self.Equilibration = EQUILIBRATION_SETTINGS
        
        self.check_valid_ff()
        self.init_parameters()
        return

    def check_valid_ff(self):

        print "\n\t\t--- AMBER DRIVER ---\n"
        print "ForceField: \t %s"    % self.forcefield
        print "Water Model:\t %s"    % self.water
        print "Mode:       \t %s\n"  % self.mode
        
        return

    def init_parameters(self):

        self.run = {} # this dict will hold the actual run parameters
        
        required_parameters = ['timestep', 'total_time', 'log_freq', 'xtc_freq',
                               'temperature', 'thermostat', 'box_size', 'Cl',
                               'Na', 'total_steps', 'barostat', 'pressure']

        # Convert time units into steps
        self.Production['total_steps'] = \
            int(self.Production['total_time'] / self.Production['timestep'])
        self.Equilibration['total_steps'] = \
            int(self.Equilibration['total_time'] / self.Equilibration['timestep'])

        # copy parameters over to 'run'
        for p in required_parameters:
            
            if self.mode == 'Equilibration': 
                if p in self.Equilibration.keys():
                    self.run[p] = self.Equilibration[p]
                else:
                    self.run[p] = self.Production[p]
                    
            elif self.mode == 'Production':
                if p == 'barostat':
                    self.run['barostat'] = 'No'
                elif p == 'pressure':
                    self.run['pressure'] = 1
                else:
                    self.run[p] = self.Production[p]
            else:
                print "Cannot understand mode: %s" % self.mode
                raise Exception("Mode must be 'Production' or 'Equilibration'")

        # check we got everything we need
        for p in required_parameters:
            assert p in self.run.keys()

        print "\nPreparing to execute %d timesteps of %s..." % (self.run['total_steps'], self.mode)

        return


def num_cores():
    result = os.sysconf('SC_NPROCESSORS_ONLN')
    return result


def check_amber_installed():

    null = open('/dev/null','w')
    r = subprocess.call("pmemd.MPI -h", shell=True, stdin=null, stdout=null, stderr=null)
    null.close()
    
    if r == 0:
        print "\nFound AMBER installation."
    else:
        raise Exception("Could not find AMBER! Looking for: pmemd.MPI in $PATH")
    
    return
    
    
def check_mpirun_installed():

    null = open('/dev/null','w')
    r = subprocess.call("mpirun -h", shell=True, stdin=null, stdout=null, stderr=null)
    null.close()

    if r == 0:
        print "\nFound MPI installation."
    else:
        raise Exception("Could not find MPI! Looking for: mpirun in $PATH")

    return


# --- Functions for writing parameter files ---

def write_min_in(Parameters):

    txt = """GENERATED BY MSMACCELERATOR: amber_driver.py
&cntrl
 imin=1,maxcyc=1000,ncyc=500,
 cut=8.0,ntb=1,
 ntc=2,ntf=2,
 ntpr=100,ntr=1
/
"""

    f = open('min.in', 'w')
    f.write(txt)
    f.close()
    print "Generated: min.in"

    return


def write_heat_in(Parameters):

    txt = """GENERATED BY MSMACCELERATOR: amber_driver.py
&cntrl
 imin=0,irest=0,ntx=1,
 nstlim=25000,dt=%f,
 ntc=2,ntf=2,
 cut=8.0, ntb=1,
 ntwx=%d,
 ntt=3, gamma_ln=2.0,
 tempi=0.0, temp0=%f,
 ntr=1,
 nmropt=1, ig=-1
/
&wt TYPE='TEMP0', istep1=0, istep2=25000,
 value1=0.1, value2=%f, /
&wt TYPE='END' /
""" % ( Parameters.run['timestep'],
        Parameters.run['xtc_freq'],
        Parameters.run['timestep'],
        Parameters.run['temperature'],
        Parameters.run['temperature'] )

    f = open('heat.in', 'w')
    f.write(txt)
    f.close()
    print "Generated: heat.in"

    return


def write_equil_in(Parameters):
    """
    Runs langevin dynamics with berendsen barostat 
    (berendsen is all thats in AMBER)
    """
    
    txt = """GENERATED BY MSMACCELERATOR: amber_driver.py
&cntrl
 imin=0,irest=0,ntx=5,
 nstlim=%d,dt=%f,
 ntc=2,ntf=2,
 cut=8.0,ntb=2,ntp=1,taup=1.0,pres0=%f
 ntwx=%d,
 ntt=3,gamma_ln=2.0,
 temp0=%f,ig=-1
/
""" % ( Parameters.run['total_steps'],
        Parameters.run['timestep'],
        Parameters.run['pressure'],
        Parameters.run['xtc_freq'],
        Parameters.run['temperature'] )

    f = open('equil.in', 'w')
    f.write(txt)
    f.close()
    print "Generated: equil.in"

    return


def write_prod_in(Parameters):
    """
    Runs langevin dynamics, NVT
    """

    txt = """GENERATED BY MSMACCELERATOR: amber_driver.py
&cntrl
 imin=0,irest=0,ntx=5,
 nstlim=%d,dt=%f,
 ntc=2,ntf=2,
 cut=8.0,ntb=1,ntp=0,
 ntwx=%d,
 ntt=3,gamma_ln=2.0,
 temp0=%f,ig=-1
/
""" % ( Parameters.run['total_steps'],
        Parameters.run['timestep'],
        Parameters.run['xtc_freq'],
        Parameters.run['temperature'] )

    f = open('equil.in', 'w')
    f.write(txt)
    f.close()
    print "Generated: equil.in"

    return


# --- Utilities

def shell_out(commands):
    """ executes commands intelligently via subprocess - waits for serial completion """
    
    log = open('driver.log','w')
    for command in commands:
        print "Executing: %s" % command
        x = subprocess.Popen( command, bufsize=-1, shell=True, stdout=log, stderr=log )
        x.wait()
        if x.returncode != 0:
            raise Exception("Error: '%s' quit with exit status: %d" % (command, x.returncode) )
    log.close()
    
    return


def clean_working_directory():

    print "Cleaning up..."
    os.system('rm \#*')
    os.system('mkdir logs')
    os.system('mv *.log logs')
    os.system('rm *.trr *.top *.itp *.edr *.cpt *.tpr out.gro conf.gro')
    os.system('mv equilibration.xtc *.mdp *.gro logs')
    
    return
    
# ----


def build_box(Parameters):



    return cmds


def run_equil():
    
    cmds = [
    "mpirun -np %d pmemd.MPI -O -i min.in -o min.out -p aln-4dua-model-complex.solv.prmtop -c aln-4dua-model-complex.solv.inpcrd -r min.rst -ref aln-4dua-model-complex.solv.inpcrd"
    "mpirun -np %d pmemd.MPI -O -i heat.in -o heat.out -p aln-4dua-model-complex.solv.prmtop -c min.rst -r heat.rst -x heat.mdcrd -ref min.rst"
    "mpirun -np %d pmemd.MPI -O -i equil.in -o equil.out -p aln-4dua-model-complex.solv.prmtop -c heat.rst -r equil.rst -x equil.mdcrd"
    ]

    return cmds


def run_md(n_threads, Parameters, wet_xtc_fn, dry_xtc_fn=None, last_wet_snapshot_fn=None):
    


    return cmds





def main(pdb_fn, forcefield, water, mode, n_threads, wet_xtc_fn, dry_xtc_fn, last_wet_snapshot_fn):

    check_amber_installed()
    check_mpirun_installed()

    print "Running on %d threads" % n_threads
    p = Parameters(forcefield, water, mode)

    # build a list of commands to execute, then shell them out
    cmds = []

    if p.mode == 'Production':
        cmds.append('pdb2gmx -f %s -ff %s -water %s' % (pdb_fn, p.forcefield, p.water))
        cmds.append('cp conf.gro start.gro')       # feed pdb2gmx out straight into grompp
        cmds.extend( run_md(n_threads, p, wet_xtc_fn, dry_xtc_fn, last_wet_snapshot_fn) )

    elif p.mode == 'Equilibration':
        # first equilibrate
        cmds.append('pdb2gmx -f %s -ff %s -water %s -ignh' % (pdb_fn, p.forcefield, p.water))
        cmds.extend( build_box(p) )
        cmds.extend( run_md(n_threads, p, 'equilibration.xtc') )

        # then seamlessly proceed to a production run
        p.mode = 'Production'
        p.init_parameters()
        cmds.append('cp end.gro start.gro')
        cmds.extend( run_md(n_threads, p, wet_xtc_fn, dry_xtc_fn, last_wet_snapshot_fn) )

    print "\nEXECUTING COMMANDS:"
    shell_out(cmds)
    clean_working_directory()

    return


if __name__ == '__main__':
    print sys.argv
    print "python %s <pdb_fn> <ff> <water> <mode> <threads>" % sys.argv[0]
    if len(sys.argv) < 6:
        raise Exception('too few arguments (sys.argv) passed')

    pdb_fn     = sys.argv[1]
    forcefield = sys.argv[2]
    water      = sys.argv[3]
    mode       = sys.argv[4]
    n_threads  = min( int(sys.argv[5]), num_cores() )

    dry_dcd_fn = 'production_wet.dcd'
    wet_dcd_fn = 'production_dry.dcd'
    last_wet_snapshot_fn = 'last_wet_snapshot.pdb'
    
    main(pdb_fn, forcefield, water, mode, n_threads,
         dry_dcd_fn, wet_dcd_fn, last_wet_snapshot_fn)


