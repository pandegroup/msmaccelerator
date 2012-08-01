import threading
import time
from Queue import Queue
from collections import namedtuple
from workqueue import WorkQueue, Task

class QMaster(threading.Thread):
    def __init__(self, project_rootdir):
        threading.Thread.__init__(self)
        
        self._job_queue = Queue()
        self.wq = WorkQueue()
        self.loop_interval = 1 # second
        self.project_rootdir = project_rootdir
    
    # Results retrieved from the workers will be saved in
    # projectroot/RawTrajectories/FORCEFIELD/job_name.xtc
    #
    # tpr files to be sent to the workers to run are saved in
    # projectrootdir/WQInputFiles/job_name.tpr
    
    def run(self):
        while True:
            if not self._job_queue.empty():
                job = self._job_queue.get()
                # save job to disk as a single input_file
                task = Task(job['name'])
                input_file = self.generate_tpr(job)
                output_file = self.output_filename(job)
                task.specify_input_file(input_file)
                task.specify_output_file(output_file)
                wq.submit(task)
            else:
                print 'Currently no jobs to submit to wq'
            
            if not wq.empty():
                task = wq.wait(self.loop_interval)
                if task:
                    # convert the output xtc to lh5
                    print 'task finished!'
            else:
                print 'no jobs in wq'
                time.sleep(self.loop_interval)
    
    def push(self, job):
        """add job to the queue"""
        print 'adding to _job_queue: %s' % job
        self._job_queue.put(job)
        
    
    def clear():
        """clear all the jobs from the queue"""
        print 'emptying queue'
        while not self.job_queue.empty():
            self._job_queue.get()
            
    def generate_tpr(self, job):
        """Create tpr file for the job and return the path to it"""
        #path = os.path.join(self.porject_rootdir, 'WQInputFiles', '%s.wqi' % job['name'])
        job_dir = os.path.join(self.project_rootdir, 'WQInputFiles', job['name'])
        os.makedirs(job_dir)
        
        job['conformation'].SaveToPDB(os.path.join(job_dir, 'structure.pdb'))
        os.copy(job['mdp_path'], os.path.join(job_dir, 'settings.mdp'))
        # run pdb2gmx to get topology and gro file
        # modify the mdp if this job contains any special parameters?
        # run grompp to get the tpr
        # save the tpr as
        job_fn = os.path.join(self.project_rootdir, 'WQInputFiles', '%s.tpr' % job['name'])
    
    def output_filename(job):
        """Return the path where the results of this job (xtc) should be put"""
        return os.path.join(self.project_rootdir, 'RawTrajectories', job['FF'], '%s.xtc' % job['Name'])
        



class Builder(object):
    def __init__(self, project_rootdir, num_states, lag_time):
        self.project_rootdir = project_rootdir
        self.num_state = num_states
        self.lag_time = lag_time
        # etc
        
    def cluster(self, round):
        # load up all of the data and cluster it, saving assignments to disk
        # Should be saved to project_rootdir/Round%d/Data/Assignments.h5
        
        # os.system(Cluster.py)
        
    def run_round(self):
        # make a new project_rootdir/Round%d directory
        # make a Trajectories folded inside that contains softlinks to 
        # the actual lh5 trajectores 
        
        # create n+1 ProjectInfo files, one of them which points to all of the
        # trajectories and then n of them, where n is the number of force fields
        # being used, that each point only to the trajectories from that force
        # field
        
        # each of these ProjectInfo files should be in a seperate directory
        # so we have:
        
        #  project_rootdir/Round2
        #           AllFF/
        #               ProjectInfo.h5 -> points to trj0.lh5, trj1.lh5, trj2.lh5
        #               Data/Assignments.h5
        #           Amber96/
        #               ProjectInfo.h5 -> points to trj0.lh5 and trj1.lh5
        #               Data/Assignments.h5 -> formed from the AllFF/Assignments.h5 by only taking the rows corresponding to the trajectories that this FF owns
        #               Data/tProb.mtx
        #               Data/tCounts.mtx
        #               etc...
        #           Amber99/
        #               ProjectInfo.h5 -> points to trj2.lh5 only
        #               Data/tProb.mtx
        #               Data/tCounts.mtx
        #               etc...
        #           Trajectories/
        #               trj0.lh5, trj1.lh5, trj2.lh5
        
        # run clustering on the AllFF ProjectInfo
        # run build_msm for each of the subdirectories
        
    
    def build_msm(self, round):
        for i in range(1, num_forcefields):
            #os.system(BuildMSM.py)
        
    
    def is_sufficient_new_data():
        # check to see if enough new lh5s are around to justify
        # building a new msm. if so return true

class SamplingBrain(object):
    def __init__(self, project_rootdir):
        self.project_rootdir = project_rootdir
    
    def generate_jobs(self):
        # load up all of the tProb.mtx from the latest round
        # pick some states to seed from based on genious statistics
        # pick some conformations from those states
        # pick a forcefield
        
        job1['Name'] = A99_S1720
        job1['FF'] = 'Amber99'
        job1['mdp'] = os.path.join(project_rootdir, '%s.mdp', ff)
        job2 = ...
        
        
        return [job1, job2]
        


def main():
    project_rootdir = os.path.abspath('.')
    qmaster = QMaster(project_rootdir)
    builder = Builder(project_rootdir, num_states, lag_time)
    brain = SamplingBrain(project_rootdir)
    
    sleep_time = 60 * 60 # 1 hour
    
    initial_round = 1
    # look in project_rootdir for the most recent Round and add 1
    
    # main loop
    for i in count(initial_round)
        jobs = brain.generate_jobs()
        for job in jobs:
            qmaster.push(job)
        
        while not builder.is_sufficient_new_data():
            time.sleep(sleep_time)
        
        builder.run_round(i)
    
    
    
