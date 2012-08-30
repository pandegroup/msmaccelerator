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

#!/usr/bin/env python

import sys
import os
import subprocess

"""
Prototype driver class, for GROMACS. This script is the part of Monakos that gets
shipped off to workers so they can actually execute MD runs. It contains all the
information needed to run any GROMACS simulation that might be required. This
should allow for the following:

 $ python gromacs_driver.py <pdb_fn> <forcefield> <water> <mode> <threads> [xtcout_fn]

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
    'thermostat'     : 'Nose-Hoover', # type
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

        print "\n\t\t--- GROMACS DRIVER ---\n"
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


def check_gromacs_installed():

    null = open('/dev/null','w')
    r = subprocess.call("grompp -h", shell=True, stdin=null, stdout=null, stderr=null)
    null.close()
    
    if r == 0:
        print "\nFound GROMACS installation."
    else:
        raise Exception("Could not find GROMACS!")
    
    return


def write_mdp(f_name, Parameters):

    txt = """; GENERATED BY MONAKOS: gromacs_driver.py
integrator               = md
dt                       = %f
nsteps                   = %d
nstxout                  = 0
nstvout                  = 0
nstlog                   = %d
nstenergy                = 0
nstxtcout                = %d
xtc_grps                 = System
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
periodic_molecules       = no
rlist                    = 1.5
coulombtype              = PME
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
optimize_fft             = yes
pme_order                = 4
fourierspacing           = 0.08
rcoulomb                 = 1.5
vdwtype                  = shift
rvdw                     = 1.25
rvdw_switch              = 1.0
tcoupl                   = %s
tc_grps                  = System
tau_t                    = 1
ref_t                    = %f
Pcoupl                   = %s
tau_p                    = 1
pcoupltype               = isotropic
compressibility          = 0.000045
ref_p                    = %f
gen_vel                  = yes
gen_temp                 = %f
constraints              = hbonds
continuation             = yes
morse                    = no
implicit_solvent         = no
""" % ( Parameters.run['timestep'],
        Parameters.run['total_steps'],
        Parameters.run['log_freq'],
        Parameters.run['xtc_freq'],
        Parameters.run['thermostat'],
        Parameters.run['temperature'],
        Parameters.run['barostat'],
        Parameters.run['pressure'],
        Parameters.run['temperature'] )

    f = open(f_name, 'w')
    f.write(txt)
    f.close()
    print "Generated: %s" % f_name

    return


def write_minimization_mdp(f_name):

    txt = """; GENERATED BY MONAKOS: gromacs_driver.py
integrator               = steep
emstep                   = 0.001
emtol                    = 10.0
nsteps                   = 100000
nstxout                  = 0
nstvout                  = 0
nstlog                   = 10
nstenergy                = 0
nstxtcout                = 0
xtc_grps                 = System
energygrps               = 
nstlist                  = 1
ns_type                  = grid
pbc                      = xyz
periodic_molecules       = no
rlist                    = 1.5
rcoulomb                 = 1.5
rvdw                     = 1.5
tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
constraints              = none
continuation             = yes
morse                    = no
implicit_solvent         = no
"""
    f = open(f_name, 'w')
    f.write(txt)
    f.close()
    print "Generated: %s" % f_name

    return


def build_box(Parameters):

    write_minimization_mdp('Minimization.mdp')

    # choose a water topology file
    w_top_dict  = {'tip3p' : 'spc216.gro'}
    if Parameters.water in w_top_dict.keys():
        water_topol = w_top_dict[Parameters.water]
    else:
        print "Could not find a topology in 'w_top_dict' for water: %s" % Parameters.water
        raise Exception("Water model has no known corresponding GROMACS topology. Please specify.")

    # format the water string
    ion_str = ''
    if Parameters.run['Cl'] > 0:
        ion_str += '-nn %d ' % Parameters.run['Cl']
    if Parameters.run['Na'] > 0:
        ion_str += '-np %d ' % Parameters.run['Na']

    cmds =[
    'editconf -f conf.gro -bt cubic -box %f %f %f -align 1 1 1' %((Parameters.run['box_size'],)*3),
    'genbox -cp out.gro -cs %s -p topol.top' % water_topol,
    'grompp -f Minimization.mdp -c out.gro -p topol.top',
    'echo SOL | genion  -s topol.tpr -o out.gro -p topol.top %s' % ion_str,
    'grompp -f Minimization.mdp -c out.gro -p topol.top',
    'mdrun -s topol.tpr -x minimization.xtc -c start.gro -pd -g EM.log'
    ]

    return cmds


def run_md(n_threads, Parameters, wet_xtc_fn, dry_xtc_fn=None, last_wet_snapshot_fn=None):
    
    write_mdp('%s_parameters.mdp' % Parameters.mode, Parameters)
    cmds = ['grompp -f %s_parameters.mdp -c start.gro -p topol.top -maxwarn 1' % Parameters.mode,
            'mdrun -nt %d -s topol.tpr -x %s -c end.gro -g prod.log' % (n_threads, wet_xtc_fn) ]

    # if running in production mode, trjconv to create a dry XTC
    if Parameters.mode == 'Production':
        assert dry_xtc_fn is not None and last_wet_snapshot_fn is not None
        cmds.append('echo 0 | trjconv -f end.gro -s topol.tpr -o %s -pbc whole' % last_wet_snapshot_fn)
        cmds.append('echo PROTEIN | trjconv -f %s -s topol.tpr -o %s -pbc whole' % (wet_xtc_fn, dry_xtc_fn))

    return cmds


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


def main(pdb_fn, forcefield, water, mode, n_threads, wet_xtc_fn, dry_xtc_fn, last_wet_snapshot_fn):

    check_gromacs_installed()

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
    if len(sys.argv) < 6: sys.exit(1)

    pdb_fn     = sys.argv[1]
    forcefield = sys.argv[2]
    water      = sys.argv[3]
    mode       = sys.argv[4]
    n_threads  = n_threads = min( int(sys.argv[5]), num_cores() )

    dry_xtc_fn = 'production_wet.xtc'
    wet_xtc_fn = 'production_dry.xtc'
    last_wet_snapshot_fn = 'last_wet_snapshot.pdb'
    
    main(pdb_fn, forcefield, water, mode, n_threads,
         dry_xtc_fn, wet_xtc_fn, last_wet_snapshot_fn)


