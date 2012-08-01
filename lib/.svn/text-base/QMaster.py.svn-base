from work_queue import Task, WorkQueue, set_debug_flag
from work_queue import WORK_QUEUE_SCHEDULE_FCFS
from work_queue import WORK_QUEUE_SCHEDULE_FILES
import time
from glob import glob
import yaml
import os, sys
import shutil
import numpy as np
import threading
import cPickle as pickle
from Queue import LifoQueue
import subprocess
import logging

from msmbuilder import Trajectory

#set_debug_flag('debug')
#set_debug_flag('wq')

class QMaster(threading.Thread):
    def __init__(self, project, port, log_freq=600): # 600 seconds
        """Initialize the qmaster
        
        Arguments:
        project_rootdir: string with the path to the project's root directory
        log_status_frequency: integer: frequency to print info about the status of the work
            queue. In units of seconds. Default is to print every 10 minutes.
        """
        threading.Thread.__init__(self)
        self.project = project
        self.log_freq = log_freq  # print time in seconds
        self.wake_freq = 1 # seconds
        
        self.wq = WorkQueue(port, name='MSMAccelerator', catalog=True, exclusive=False)
        self.logger = logging.getLogger('MSMAccelerator.QMaster')
        self.logger.info('WORK QUEUE MASTER LISTENING ON PORT: {0}'.format(self.wq.port))
        self.logger.info('(Start a local worker with >> work_queue_worker -d all localhost {0} & )'.format(self.wq.port))
        
        # method controls whether or not we need to bring back solvated_xtc as well
        if self.project.method == 'explicit':
            self.return_wet_xtc = True
        elif self.project.method == 'implicit':
            self.return_wet_xtc = False
        else:
            raise 
        self.logger.info('Return wet xtc set to %s', self.return_wet_xtc)
        
        # what does this specify algorithm do?
        self.wq.specify_algorithm(WORK_QUEUE_SCHEDULE_FCFS)
        
        # fast abort kills jobs that appear to be stragling (taking more than 1.5x average)
        #self.wq.activate_fast_abort(1.5)
        
        # setting the stop event signals for the thread to die
        self._stop = threading.Event()
        
        # the thread sets the  event every time a job returns or there are no waiting jobs
        # and it finished post processing. See the wait method
        self._mainloop_wake_event_cause = None
        self._mainloop_wake_event = threading.Event()
        
        # start the thread
        self.start()
    
    def run(self):
        """Main thread-loop for the QMaster thread"""
        last_print = time.time()
        
        while True:
            time.sleep(self.wake_freq)
            
            if not self.wq.empty():
                t = self.wq.wait(self.wake_freq)
                if t:
                    job_dict = pickle.loads(t.tag)
                    job_dict['returned_time'] = time.time()
                    job_dict['host'] = t.host
                    if t.return_status != 0:
                        self.logger.error('Worker returned nonzero exit status for job: {0}'.format(str(job_dict)))
                    else:
                        self.on_return(job_dict)
                    self._mainloop_wake_event_cause = 'job returned'
                    self._mainloop_wake_event.set()
            
            if self.wq.stats.tasks_waiting == 0 and not self._mainloop_wake_event.is_set():
                self._mainloop_wake_event_cause = 'queue empty'
                self._mainloop_wake_event.set() # also set the event if there are no tasks in the queue

            if self._stop.is_set():
                self.logger.info('Recieved stop signal. Shutting down all workers')
                self.wq.shutdown_workers(0) # 0 indicates to shut all of them down
                sys.exit(0)
            
            if time.time() - last_print > self.log_freq:
                self.logger.info('workers initialized: %d, ready: %d, busy: %d' % (self.wq.stats.workers_init, self.wq.stats.workers_ready, self.wq.stats.workers_busy))
                self.logger.info('workers running: %d, waiting: %d, complete: %d' % (self.wq.stats.tasks_running, self.wq.stats.tasks_waiting, self.wq.stats.tasks_complete))
                last_print = time.time()

    def num_jobs_waiting(self):
        return self.wq.stats.tasks_waiting

    def num_jobs_in_queue(self):
        """Get the number of jobs currently in the work queue"""
        return self.wq.stats.tasks_running + self.wq.stats.tasks_waiting
        
    def stop(self):
        """Signal the Qmaster thread to stop"""
        self._stop.set()

    def wait(self):
        self._mainloop_wake_event.wait()
        self._mainloop_wake_event.clear()
        cause = self._mainloop_wake_event_cause
        assert cause in ['job returned', 'queue empty']
        return cause
    
    def submit(self, job):
        """ Submit a job to the work-queue for further sampling.
        
        Job should be a dict containing at least the following fields
        
        'name': The name of the job
        'conf': MSMBuilder Conformation object with the starting structure
        'ff': Name of the forcefield to use. (string)
        'water': Name of the water model to use (string)
        
        'mode': 'Production' or 'Equilibration'
        'threads': Number of threads to use on target machine (int)
        'driver': path to MD driver script (e.g. Monakos' gromacs_driver.py)
        'output_extension': path to 
        """
        # if job is a list, submit each individually
        if isinstance(job, list):
            # submit the new jobs. We're submitting them to a LIFO queue,
            # so we submit them in reverse order so that the first one
            # the client submitted is on top
            for j in job[::-1]:
                self.submit(j)
            return
        
        # make sure that job contains the right fields
        required = ['name', 'conf', 'ff', 'water', 'mode', 'threads', 'driver']
        def confirm(item):
            if not item in job:
                raise KeyError('job needs to contain the key "%s"' % item)
        map(confirm, required)
        
        fileroot = str(np.random.randint(sys.maxint))
        job['fileroot'] = fileroot

        pdb_fn =  '%s.pdb' % fileroot
        dry_xtc_fn = '%s_dry.xtc' % fileroot
        wqlog_fn = '%s.wqlog' % fileroot
        
        xtc_dir = self.project.xtc_dir(job['ff'])
        job['dry_xtc'] = os.path.join(xtc_dir, dry_xtc_fn)
        job['last_wet_snapshot'] = os.path.join(self.project.last_snapshot_dir(job['ff']), pdb_fn)
        job['wqlog'] = os.path.join(self.project.wqlog_dir(job['ff']), wqlog_fn)
        job['pdb_fn'] = os.path.join(self.project.starting_conf_dir(job['ff']), pdb_fn)
        shutil.move(job['conf'], job['pdb_fn'])
        del job['conf']
        
        job['submit_time'] = time.time()
        
        driver_fn = os.path.split(job['driver'])[1]
        
        task = Task('python ./%s %s %s %s %s %s > %s' % \
                    (driver_fn, pdb_fn, job['ff'], job['water'], \
                     job['mode'], job['threads'], wqlog_fn))
        
        task.specify_input_file(job['driver'], driver_fn)
        task.specify_input_file(job['pdb_fn'], pdb_fn)
        task.specify_output_file(job['wqlog'], 'logs/driver.log')
        
        # this is the xtc without water
        remote_output_fn = 'production_dry%s' % job['output_extension']
        task.specify_output_file(job['dry_xtc'], remote_output_fn)
        
        if self.return_wet_xtc:
            wet_xtc_fn = '%s_wet.xtc' % fileroot
            # this is the XTC file with waters, generated by the driver
            # when you're doing implicit solvent only, this stuff is not used.
            remote_output_fn = 'production_wet%s' % job['output_extension']
            job['wet_xtc'] = os.path.join(xtc_dir, wet_xtc_fn)
            task.specify_output_file(job['wet_xtc'], remote_output_fn)
            task.specify_output_file(job['last_wet_snapshot'], 'last_wet_snapshot.pdb')
        else:
            self.logger.debug('Not requesting production_wet%s from driver (implicit)' % job['output_extension'])
        
        # pickle the job dict and set it as the tag of the wq task
        task.specify_tag(pickle.dumps(job))
        task.specify_algorithm(WORK_QUEUE_SCHEDULE_FILES) # what does this do?
        self.wq.submit(task)
        
        self.logger.info('Submitted to queue: {0} ({1}, {2})'.format(fileroot, job['ff'], job['name']))
        
    def on_return(self, job):
        """Called by main thread on the return of data from the workers.
        Post-processing"""
        self.logger.info('Retrieved "{0}" xtc. Converting to lh5...'.format(job['name']))
        
        try:
            # save lh5 version of the trajectory
            traj_dir = self.project.traj_dir(job['ff'])
            trajnum = len(glob(os.path.join(traj_dir, '*.lh5')))
            lh5_fn = os.path.abspath(os.path.join(traj_dir, '%d.lh5' % trajnum))
            conf = Trajectory.LoadTrajectoryFile(self.project.pdb_topology_file)
            traj = Trajectory.LoadTrajectoryFile(job['dry_xtc'], Conf=conf)
            traj.SaveToLHDF(lh5_fn)
        
        except Exception as e:
            self.logger.error('When postprocessing {0}, convert to lh5 failed!'.format(str(job)))
            self.logger.exception(e)
            raise
        
        # create softlink to the lh5 trajectory in the JointFF directory
        softlink_dir = self.project.traj_dir(self.project.joint_ff['name'])
        
        softlink_num = len(glob(os.path.join(softlink_dir, '*.lh5')))
        softlink_fn = os.path.join(softlink_dir, '%d.lh5' % softlink_num)
        os.symlink(lh5_fn, softlink_fn)

        # update the TrajLog file
        job['AllFF_fn'] = softlink_fn
        job['lh5_fn'] = lh5_fn
        job['TrajLength'] = len(traj)
        job['lh5_trajnum'] = trajnum
        self.project.traj_log_add(job)
        self.logger.info('Finished converting new traj to lh5 sucessfully')
    
if __name__ == '__main__':
    q = QMaster('.')
    job = {'name': 'testjob3',
           'driver': 'python /home/rmcgibbo/monakos/drivers/GROMACS/gromacs_driver.py',
           'conf': Trajectory.LoadFromPDB('/home/rmcgibbo/monakos/drivers/GROMACS/ala5.pdb'),
           'ff': 'amber99sb-ildn',
           'water': 'tip3p',
           'mode': 'equilibration',
           'threads': 8}
    
    q.submit(job)
    q.stop()
