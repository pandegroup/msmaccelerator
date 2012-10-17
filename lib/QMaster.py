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

from work_queue import Task, WorkQueue, set_debug_flag
from work_queue import WORK_QUEUE_SCHEDULE_FCFS
from work_queue import WORK_QUEUE_SCHEDULE_FILES
import time
import os, sys
import threading
import cPickle as pickle
import logging
from datetime import datetime

import msmbuilder.Trajectory
import models
from database import Session, with_db_lock
from utils import save_file, load_file


logger = logging.getLogger('MSMAccelerator.QMaster')
#set_debug_flag('debug')
#set_debug_flag('wq')

class QMaster(threading.Thread):
    def __init__(self, project, port, log_freq=600): # 600 seconds
        """Initialize the QMaster
        
        Parameters
        ----------
        project : 
        port : int
        log_freq : int, optional
            frequency to print info about the status of the work queue.
            In units of seconds. Default is to print every 10 minutes.
        """
        
        threading.Thread.__init__(self)
        self.project = project
        
        self.log_freq = log_freq  # print time in seconds
        self.wake_freq = 1 # seconds
        
        self.wq = WorkQueue(port, name='MSMAccelerator', catalog=True, exclusive=False)

        logger.info('WORK QUEUE MASTER LISTENING ON PORT: %d', self.wq.port)
        logger.info('(Start a local worker with >> work_queue_worker -d all localhost %d & )', self.wq.port)
        
        # method controls whether or not we need to bring back solvated_xtc as well
        if self.project.method == 'explicit':
            self.return_wet_xtc = True
        elif self.project.method == 'implicit':
            self.return_wet_xtc = False
        else:
            raise Exception("project.method must be 'explicit' or 'implicit'")
        logger.info('Return wet xtc set to %s', self.return_wet_xtc)
        
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
                    if t.return_status != 0:
                        logger.error('Worker returned nonzero exit status for job: %d', t.return_status)
                    else:
                        self.on_return(t)
                    self._mainloop_wake_event_cause = 'job returned'
                    self._mainloop_wake_event.set()
            
            if self.wq.stats.tasks_waiting == 0 and not self._mainloop_wake_event.is_set():
                self._mainloop_wake_event_cause = 'queue empty'
                self._mainloop_wake_event.set() # also set the event if there are no tasks in the queue

            if self._stop.is_set():
                logger.info('Recieved stop signal. Shutting down all workers')
                self.wq.shutdown_workers(0) # 0 indicates to shut all of them down
                sys.exit(0)
            
            if time.time() - last_print > self.log_freq:
                logger.info('workers initialized: %d, ready: %d, busy: %d', self.wq.stats.workers_init, self.wq.stats.workers_ready, self.wq.stats.workers_busy)
                logger.info('workers running: %d, waiting: %d, complete: %d', self.wq.stats.tasks_running, self.wq.stats.tasks_waiting, self.wq.stats.tasks_complete)
                last_print = time.time()

    def num_jobs_waiting(self):
        """Number of jobs waiting to be sent out
        
        This number should be kept at 1, and when it drops to zero a new job
        should be generated.
        
        Returns
        -------
        n : int
            The number
        """
        return self.wq.stats.tasks_waiting

    def num_jobs_in_queue(self):
        """Get the number of jobs currently in the work queue
        
        This includes both the jobs running remotely and the ones waiting
        here
        
        Returns
        -------
        n : int
            The number
        """
        return self.wq.stats.tasks_running + self.wq.stats.tasks_waiting
        
    def stop(self):
        """Signal the Qmaster thread to stop"""
        self._stop.set()

    def wait(self):
        """Block until some sort of action happens in the main-thread loop.
        
        This call will return either when a job as returned from the workers,
        or when the queue is empty (last job in the local queue has been sent
        out)
        
        Returns
        -------
        wakeup_cause : str
            Either 'job returned' or 'queue empty', depending on the reason
        """
        self._mainloop_wake_event.wait()
        self._mainloop_wake_event.clear()
        cause = self._mainloop_wake_event_cause
        if not cause in ['job returned', 'queue empty']:
            raise Exception('Bad wakeup cause')
        return cause
    
    @with_db_lock
    def submit(self, traj):
        """ Submit a job to the work-queue for further sampling.
        
        Parameters
        ----------
        """
        if traj.submit_time is not None:
            raise ValueError("This traj has already been submitted")
        Session.add(traj)
        Session.flush()
        traj.populate_default_filenames()
        
        if not hasattr(traj, 'init_pdb'):
            raise ValueError('Traj is supposed to have a pdb object tacked on')            
        save_file(traj.init_pdb_fn, traj.init_pdb)
        
        remote_driver_fn = os.path.split(str(traj.forcefield.driver))[1]
        remote_pdb_fn = 'input.pdb'
        remote_output_fn = 'production_dry{}'.format(traj.forcefield.output_extension)
        
        if traj.mode is None or traj.forcefield is None:
            raise ValueError('malformed traj')

        task = Task('python ./{driver} {pdb_fn} {ff} {water} {mode} {threads}'.format(
            pdb_fn=remote_pdb_fn,
            mode=traj.mode,
            driver=remote_driver_fn,
            ff=traj.forcefield.name,
            water=traj.forcefield.water,
            threads=traj.forcefield.threads))
        
        
        #why does traj.forcefield.driver come out as unicode?
        task.specify_input_file(str(traj.forcefield.driver), remote_driver_fn)
        task.specify_output_file(traj.wqlog_fn, 'logs/driver.log')
        task.specify_input_file(traj.init_pdb_fn, remote_pdb_fn)
        task.specify_output_file(traj.dry_xtc_fn, remote_output_fn)
        
        if self.return_wet_xtc:
            # this is the XTC file with waters, generated by the driver
            # when you're doing implicit solvent only, this stuff is not used.
            remote_wet_output_fn = 'production_wet{}'.format(traj.forcefield.output_extension)
            task.specify_output_file(traj.wet_xtc_fn, remote_wet_output_fn)
            task.specify_output_file(traj.last_wet_snapshot_fn, 'last_wet_snapshot.pdb')
        else:
            logger.debug('Not requesting production_wet%s from driver (implicit)', traj.forcefield.output_extension)
        
        task.specify_tag(str(traj.id))
        task.specify_algorithm(WORK_QUEUE_SCHEDULE_FILES) # what does this do?
        
        traj.submit_time = datetime.now()

        self.wq.submit(task)    
        logger.info('Submitted to queue: %s', traj)
    
    @with_db_lock
    def on_return(self, task):
        """Called by main thread on the return of data from the workers.
        Post-processing"""
        logger.info('Retrieved task %s', task.tag)
        traj = Session.query(models.Trajectory).get(int(task.tag))
        
        try:
            # save lh5 version of the trajectory
            conf = load_file(self.project.pdb_topology_file)
            coordinates = msmbuilder.Trajectory.load_trajectory_file(str(traj.dry_xtc_fn), Conf=conf)
            save_file(traj.lh5_fn, coordinates)
        
        except Exception as e:
            logger.error('When postprocessing %s, convert to lh5 failed!', traj)
            logger.exception(e)
            raise
        
        # convert last_wet_snapshot to lh5
        pdb_to_lh5(traj, 'last_wet_snapshot_fn')
        pdb_to_lh5(traj, 'init_pdb_fn')


        traj.host = task.host
        traj.returned_time = datetime.now()
        traj.length = len(coordinates)
        logger.info('Finished converting new traj to lh5 sucessfully')


def pdb_to_lh5(traj, field):
    path = getattr(traj, field)
    data = load_file(path)
    new_fn = os.path.splitext(path)[0] + '.lh5'
    save_file(new_fn, data)
    os.unlink(path)
    setattr(traj, field, new_fn)
    
