
"""
This file contains all of the functionality for running MSMAccelerator in
'batch' mode.
"""

import os, sys
import time
import subprocess
from subprocess import PIPE

import shutil
import logging
from glob import glob

from msmbuilder import Trajectory

# set logging
logger = logging.getLogger('MSMAccelerator.BatchMaster')


def get_working_directory(working_dir_base):
    """ Ensures that a clean working directory is generated an available """

    # choose a directory to work in
    r = os.popen('whoami')
    username = r.read().strip()

    working_dir = os.path.join(working_dir_base, username, 'msmaccelerator-batch')

    if os.path.exists(working_dir):
        logger.info("Working in: %s" % working_dir)
    else:
        try:
            os.mkdir(working_dir)
        except:
            logger.critical("Check permissions on: %s" % working_dir_base)
            raise Exception("Cannot generate directory: %s for scratch!" % working_dir)

    return working_dir


def write_production_job_script(round_num, num_cpus, brain, job_script_name='msma-batch-production.sh', 
                                working_dir_base='/scratch', scheduler='PBS'):
    """ Generates a scheduler script that can be q-subbed to start a bunch of MD runs.
        This function first calls the brain and tells it to generate a bunch of jobs,
        and then formats a script that will start all of those jobs.

        Note that specifics about what resources to use should be specified when this
        scheduler script is actually submitted to the queue (see write_builder_job_script())
    """

    working_dir = get_working_directory(working_dir_base)

    used_cpus = 0
    num_jobs = 0
    jobs = []

    job_string = ''
    while used_cpus < num_cpus:

        job_dict = brain.generate_job(round_num, tmp_dir=working_dir)
        job_dict['fileroot'] = os.path.basename(job_dict['conf'])[:-4] # number identifier 

        job_string += 'mkdir -p {wd}/{fileroot}; cd {wd}/{fileroot}; python {driver} {conf} {ff} {water} {mode} {threads} & \n'.format(wd=working_dir, **job_dict)
        used_cpus += job_dict['threads']
        num_jobs += 1

        jobs.append( job_dict )

    logger.info("Generating %d jobs..." % num_jobs)


    if scheduler == 'PBS':
        txt = """#!/bin/bash

#PBS -N MSMA-batch-production
#PBS -e MSMA-batch-production.log
#PBS -o MSMA-batch-production.log
#PBS -V

cd {workdir}
echo working in: `pwd`

# execute commands of interest
{job_string}

wait""".format(workdir=working_dir, job_string=job_string)

    else:
        raise Exception('Scheduler option: %s not valid/implemented' % scheduler)

    f = open(job_script_name, 'w')
    f.write(txt)
    f.close()
    logger.info("Generated: %s" % job_script_name)

    return jobs


def write_builder_job_script(job_script_name='msma-batch-builder.sh', scheduler='PBS'):
    """ Writes a scheduler script that starts a builder job.
    """ 

    working_dir = get_working_directory(working_dir_base)

    job_string = "ForceBuilder %s" % project.params_fn

    # write the scheduler script
    if scheduler == 'PBS':
        txt = """#!/bin/bash

#PBS -N MSMA-batch-production
#PBS -e MSMA-batch-production.log
#PBS -o MSMA-batch-production.log
#PBS -V

cd {workdir}
echo working in: `pwd`

# execute commands of interest
{job_string}

wait""".format(workdir=working_dir, job_string=job_string)

    else:
        raise Exception('Scheduler option: %s not valid/implemented' % scheduler)

    f = open(job_script_name, 'w')
    f.write(txt)
    f.close()
    logger.info("Generated: %s" % job_script_name)

    return


def submit_job_to_queue(job_script_name, qsub_arguments, scheduler='PBS'):
    """ Simply shells out a qsub command, with some user-specified dressing """

    if scheduler == 'PBS':
        cmd = "qsub %s %s" % (qsub_arguments, job_script_name)
        logger.info("Executing: %s" % cmd)
        p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)


    else:
        raise Exception('Scheduler option: %s not valid/implemented' % scheduler)


    return



def wait_for_job_completion(sleep_time=600, scheduler='PBS'):
    """ This function checks if MSMAccelerator is running, and if so 
        waits 'sleep_time' 
    """

    if scheduler == 'PBS':

        job_running = True
        while job_running:

            time.sleep(sleep_time)

            cmd = "qstat -u `whoami` | grep MSMA-batch | wc -l"
            p = subprocess.Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
            stdout, stderr = p.communicate()
            if stderr == '':
                active_jobs = int(stdout.strip())
                logger.info("Found %d running jobs matching 'MSMA-batch-*'" % active_jobs)
            else:
                raise Exception("Error running qstat: %s %s" % (stdout, stderr))

            if active_jobs != 0:
                logger.info("Sleeping for %d seconds while we wait for it to finish" % sleep_time)
                job_running = True
            else:
                logger.info("All remote jobs appear to be finished.")
                job_running = False

    else:
        raise Exception('Scheduler option: %s not valid/implemented' % scheduler)

    return


def process_completed_trajectories(project, jobs, working_dir_base='/scratch'):
    """ Takes all of the finished jobs and processes them
    """

    logger.info("Post processing %s jobs" % len(jobs))
    working_dir = get_working_directory(working_dir_base)

    for job in jobs:

        # first, just set a bunch of metadata for the job
        fileroot = job['fileroot']

        pdb_fn =  '%s.pdb' % fileroot
        dry_xtc_fn = '%s_dry.xtc' % fileroot
        wqlog_fn = '%s.wqlog' % fileroot

        xtc_dir = project.xtc_dir(job['ff'])

        job['dry_xtc'] = os.path.join(xtc_dir, dry_xtc_fn)
        job['last_wet_snapshot'] = os.path.join( project.last_snapshot_dir(job['ff']), pdb_fn)
        job['wqlog'] = os.path.join( project.wqlog_dir(job['ff']), wqlog_fn)
        job['pdb_fn'] = os.path.join( project.starting_conf_dir(job['ff']), pdb_fn)
        shutil.move(job['conf'], job['pdb_fn'])
        del job['conf']

        driver_fn = os.path.split(job['driver'])[1]

        # next, move all of the files from the working_dir to the project directory
        job_working_dir = os.path.join( working_dir, fileroot )
        
        # copy over the driver log
        driver_log = os.path.join( job_working_dir, 'logs/driver.log')
        shutil.copy( driver_log, job['wqlog'] )
        #task.specify_output_file(job['wqlog'], 'logs/driver.log')
 
        # this is the xtc without water
        remote_output_fn = 'production_dry%s' % job['output_extension']
        shutil.copy( os.path.join( job_working_dir, remote_output_fn ), job['dry_xtc'])
        #task.specify_output_file(job['dry_xtc'], remote_output_fn)

        if project.method == 'explicit':
            wet_xtc_fn = '%s_wet.xtc' % fileroot
            # this is the XTC file with waters, generated by the driver
            # when you're doing implicit solvent only, this stuff is not used.
            remote_output_fn = 'production_wet%s' % job['output_extension']
            job['wet_xtc'] = os.path.join(xtc_dir, wet_xtc_fn)
            shutil.copy( os.path.join( job_working_dir, remote_output_fn ), job['wet_xtc'])
            shutil.copy( os.path.join( job_working_dir, 'last_wet_snapshot.pdb'), job['last_wet_snapshot'])
            #task.specify_output_file(job['wet_xtc'], remote_output_fn)
            #task.specify_output_file(job['last_wet_snapshot'], 'last_wet_snapshot.pdb')
        else:
            logger.debug('Not requesting production_wet%s from driver (implicit)' % job['output_extension'])

        # then, update the internal project directory structure
        # and update the TrajLog.yaml file
        try:
            # save lh5 version of the trajectory
            traj_dir = project.traj_dir(job['ff'])
            trajnum = len(glob(os.path.join(traj_dir, '*.lh5')))
            lh5_fn = os.path.abspath(os.path.join(traj_dir, '%d.lh5' % trajnum))
            conf = Trajectory.LoadTrajectoryFile( project.pdb_topology_file)
            traj = Trajectory.LoadTrajectoryFile(job['dry_xtc'], Conf=conf)
            traj.SaveToLHDF(lh5_fn)

        except Exception as e:
            logger.error('When postprocessing {0}, convert to lh5 failed!'.format(str(job)))
            logger.exception(e)
            raise

        # create softlink to the lh5 trajectory in the AllFF directory
        softlink_dir = project.traj_dir('AllFF')
        softlink_num = len(glob(os.path.join(softlink_dir, '*.lh5')))
        softlink_fn = os.path.join(softlink_dir, '%d.lh5' % softlink_num)
        os.symlink(lh5_fn, softlink_fn)

        # update the TrajLog file
        job['AllFF_fn'] = softlink_fn
        job['lh5_fn'] = lh5_fn
        job['TrajLength'] = len(traj)
        job['lh5_trajnum'] = trajnum
        project.traj_log_add(job)
        logger.info('Finished converting new traj to lh5 sucessfully')

        # finally, if all of that worked (whew), we want to clean up the scratch space
        logger.info('Removing directory: %s' % job_working_dir)
        shutil.rmtree( job_working_dir )

    return


