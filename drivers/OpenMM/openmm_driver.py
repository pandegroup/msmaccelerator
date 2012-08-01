
import sys
import os

from simtk.openmm.app import * 
from simtk.openmm import * 
from simtk.unit import *


# + ------------    Manually set parameters   ------------ +
#
# this constitutes a minimalist set of parameters you
# may want to use in an MD run - further control
# can be gained by directly modifying the .mdp text in the
# write_mdp method below

PRODUCTION_SETTINGS = {
    'timestep'       : 0.002,         # ps
    'total_time'     : 1,             # ps
    'dcd_freq'       : 1,           # timesteps
    'temperature'    : 300,           # Kelvin
    'box_size'       : 4,             # Angstroms
    'Cl'             : 0,             # number Cl's
    'Na'             : 0              # number Na's
}

# the equilibration run copys over any parameters not
# specified from the Production set. So just add any-
# thing that should be changed for production

EQUILIBRATION_SETTINGS = {
    'timestep'       : 0.002,         # ps
    'total_time'     : 1,             # ps
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

        print "\n\t\t--- OpenMM DRIVER ---\n"
        print "ForceField: \t %s"    % self.forcefield
        print "Water Model:\t %s"    % self.water
        print "Mode:       \t %s\n"  % self.mode
        
        return

    def init_parameters(self):

        self.run = {} # this dict will hold the actual run parameters
        
        required_parameters = ['timestep', 'total_time', 'dcd_freq',
                               'temperature', 'box_size', 'Cl',
                               'Na', 'total_steps', 'pressure']

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


def explicit_run(p, pdb, forcefield, dry_dcd_fn, wet_dcd_fn, last_wet_snapshot):

    print "Running in explicit solvent"
    integrator = VerletIntegrator( p.run['timestep'] * picoseconds )
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, \
                                     nonbondedCutoff=1.2*nanometer, constraints=HBonds)

    # if we are equilibrating, run in NPT
    if p.mode == 'Equilibration':

        print "Running Equilibration..."

        # add pressure and temperature coupling
        system.addForce( AndersenThermostat( p.run['temperature']*kelvin, 1/picosecond ) )
        system.addForce( MonteCarloBarostat( p.run['pressure']*bar, p.run['temperature']*kelvin) )
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        print "\tPlatform: %s" % (simulation.context.getPlatform().getName())

        print "\tMinimizing Energy"
        simulation.minimizeEnergy()
       
        # at the end, write a PDB for the equilibrated box to load in for production
        print "\tRunning MD"
        simulation.step(p.run['total_steps'])
        equilibrated_positions = simulation.context.getPositions()
        pdb.writeFile( simulation.topology, equilibrated_positions,
                       file='equilibrated.pdb' ) # write the equilibrated file

        # rebuild the simulation for production
        print "\tRebuilding simulation for production"
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, \
                                         nonbondedCutoff=1.2*nanometer, constraints=HBonds)

        # switch to Production mode
        p.mode = 'Production'
        integrator = VerletIntegrator( p.run['timestep'] * picoseconds ) 

    elif p.mode == 'Production':
        pass

    else:
        raise Exception('Invaid mode: %s' % p.mode)

    # then do a production run, regardless of what mode we're in
    print "Running Production..."

    system.addForce( AndersenThermostat( p.run['temperature']*kelvin, 1/picosecond ) )
    simulation.context.setPositions( equilibrated_positions )
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.reporters.append( DCDReporter(wet_dcd_fn, p.run['dcd_freq']) )
    simulation.step( p.run['total_steps'] )
    pdb.writeFile( simulation.topology, equilibrated_positions, file=last_wet_snapshot ) # write the final wet file

    # TJL to self: still need to report the dry_dcd filename

    return


def implicit_run(p, pdb, forcefield, dry_dcd_fn):

    print "Preparing implicit solvent run..."

    integrator = LangevinIntegrator( p.run['temperature']*kelvin, 91./picosecond, p.run['timestep'] * picoseconds )
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic, \
                                     nonbondedCutoff=1.0*nanometer, constraints=HBonds)

    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    print "\tPlatform: %s" % (simulation.context.getPlatform().getName())

    print "\tMinimizing Energy"
    simulation.minimizeEnergy()

    print "\tRunning MD"
    simulation.reporters.append( DCDReporter(dry_dcd_fn, p.run['dcd_freq']) )
    simulation.step( p.run['total_steps'] )

    return


def main(pdb_fn, forcefield, water, mode, n_threads, dry_dcd_fn, wet_dcd_fn, last_wet_snapshot):

    print "Running on %d threads" % n_threads

    # get the value of the $CUDA_DEVICE env. var. and use it to tell us what GPU to run on
    try:
        deviceid = os.environ['CUDA_DEVICE']
    except:
        deviceid = 0
    platform = openmm.Platform_getPlatform(1)
    platform.setPropertyDefaultValue('CudaDevice', '%d' % deviceid)
    print "Device ID set to: %d" % deviceid

    p = Parameters(forcefield, water, mode)
   
    # set up the system coordinates and topology
    print "\nLoading:\n\t%s\n\t%s\n\t%s\n" % ( pdb_fn, p.forcefield+'.xml', p.water+'.xml' )
    pdb = PDBFile( pdb_fn )
    forcefield = ForceField(p.forcefield+'.xml', p.water+'.xml')

    if p.water in ['amber03_gbvi', 'amber10_gbvi', 'amber96_gbvi', 'amber99_gbvi', 'amoeba2009_gk',
                   'amber03_obc', 'amber10_obc', 'amber96_obc', 'amber99_obc' ]:
        implicit_run(p, pdb, forcefield, dry_dcd_fn)
    else:
        explicit_run(p, pdb, forcefield, dry_dcd_fn, wet_dcd_fn, last_wet_snapshot)


    print "Driver finished gracefully"
    return


if __name__ == '__main__':
    print sys.argv
    print "python %s <pdb_fn> <ff> <water> <mode> <threads>" % sys.argv[0] 
    if len(sys.argv) < 6: sys.exit(1)

    pdb_fn     = sys.argv[1]
    forcefield = sys.argv[2]
    water      = sys.argv[3]
    mode       = sys.argv[4]
    n_threads  = min( int(sys.argv[5]), num_cores() )

    dry_dcd_fn = 'production_dry.dcd'
    wet_dcd_fn = 'production_wet.dcd'
    last_wet_snapshot = 'last_wet_snapshot.pdb'
    
    main( pdb_fn, forcefield, water, mode, n_threads,
          dry_dcd_fn, wet_dcd_fn, last_wet_snapshot )



