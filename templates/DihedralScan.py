#!/usr/bin/env python

from Settings import AllDict, OptDict, SPDict, TheoryDict, BasisDict
from leeping import molecule, qcparse
import re
import os
import itertools
import subprocess as sp
import numpy as n
import shutil
import work_queue

OBPath = "/home/leeping/opt/openbabel-2.3.1/bin"
cwd = os.getcwd()

work_queue.set_debug_flag('all')
wq = work_queue.WorkQueue(port=5813, exclusive=False, shutdown=True)
wq.specify_name('test')

def qc_section(qcopts,section):
    if section not in qcopts:
        return []
    elif qcopts[section] == "":
        return []
    else:
        Answer = ["$%s" % section]
        for line in qcopts[section]:
            Answer.append(line.strip())
        Answer.append("$end")
    return Answer

def qc_rem_section(qcopts):
    Answer = ["$rem"]
    rems = qcopts['rem']
    if 'jobtype' in rems:
        Answer.append("%-20s %s" % ('jobtype',rems['jobtype']))
    if 'exchange' in rems:
        Answer.append("%-20s %s" % ('exchange',rems['exchange']))
    if 'correlation' in rems:
        Answer.append("%-20s %s" % ('correlation',rems['correlation']))
    if 'basis' in rems:
        Answer.append("%-20s %s" % ('basis',rems['basis']))
    for option in sorted([i for i in rems]):
        if option not in ['jobtype','exchange','correlation','basis']:
            Answer.append ("%-20s %s" % (option,rems[option]))
    Answer.append("$end")
    return Answer
        
class ScanJob:
    def __init__(self,line):
        s = line.split()
        self.Number              = int(s[0])
        self.Name                = s[1]
        self.Dir                 = "%i_%s" % (self.Number,self.Name)
        self.Smiles              = s[2]
        # First element in the list is the ATOMS
        # Second element in the list are the SCAN VALUES
        self.Scans               = [(i.split('[')[0].split(','), \
                                     [int(j) for j in i.split('[')[1].split(']')[0].split(':')]) for i in s[3].split(';')]
        self.Calcs               = []
        # Set the default spacing which depends on the number of scanning dimensions
        if len(self.Scans) == 1:
            self.Spac = 3
        elif len(self.Scans) == 2:
            self.Spac = 10
        elif len(self.Scans) == 3:
            self.Spac = 15
        else:
            pause("wtf? We need to have 1-3 scan dimensions.")

    def add_calc(self,line):
        self.Calcs.append({"Charge"    : 0,
                           "Mult"      : 1,
                           "SPTheory"  : line.strip().split('//')[0].split('/')[0].upper(),
                           "SPBas"     : line.strip().split('//')[0].split('/')[1].upper(),
                           "OptTheory" : line.strip().split('//')[1].split('/')[0].upper(),
                           "OptBas"    : line.strip().split('//')[1].split('/')[1].upper()
                           })
    def scan(self):
        absdir = os.path.join(cwd,'Molecules',self.Dir)
        print "Now working on the dihedral scan for %s" % self.Name
        # Make the job directory
        mkdir(absdir)
        # Create the initial .xyz file
        self.initxyz = os.path.join(cwd,'Molecules',self.Dir,'initial.xyz')
        if not os.path.exists(self.initxyz):
            with open(os.path.join(cwd,'temp','temp.smi'),'w') as f: f.write(self.Smiles+"\n")
            sp.call('%s/babel -ismi %s/temp/temp.smi -oxyz %s/initial.xyz --title %s --gen3D' % (OBPath,cwd,absdir,self.Name),shell=True)
        self.all_vals = self.prepare_angles()
        touch(os.path.join(absdir,'scan_start.xyz'))
        for calc in self.Calcs:
            Answer = []
            AnsVec = []
            PrintOut = []
            for val in self.all_vals:
                Energy = self.evaluate(calc, val)
                if Energy != None:
                    AnsVec.append(Energy)
                Answer.append((val,Energy))
            if len(AnsVec) > 0:
                E0 = min(AnsVec)
                Fac = 627.51
                for i in Answer:
                    if i[1] != None:
                        PrintOut.append(' '.join(["%12.6f" % j for j in i[0]]) + ' % 12.6f' % ((i[1]-E0)*Fac) + '\n')
                # Note that this file name has OptTheory and OptBas come second, whereas the calculation directories have it come first
                with open(os.path.join(cwd,'Data',self.Dir+'_%s_%s__%s_%s.txt' % (calc['SPTheory'],calc['SPBas'],calc['OptTheory'],calc['OptBas'])),'w') as f: f.writelines(PrintOut)
        
    def prepare_angles(self):
        ranges = []
        for i in self.Scans:
            lo = i[1][0]
            hi = i[1][-1]
            spac = i[1][1] if len(i[1]) == 3 else self.Spac
            ranges.append(list(n.linspace(lo, hi, 1+(hi-lo)/spac)))
        return list(itertools.product(*ranges))

    def create_qcin(self, M, calc, val, outfnm, jobtype):
        ## Default general options - basically a collapsed veresion of gen_opts_types.
        qcopts = {"rem":{},"basis":""}
        RemDict = OptDict if jobtype.upper() == 'OPT' else SPDict
        LevelTheory = calc["OptTheory"] if jobtype.upper() == 'OPT' else calc["SPTheory"]
        BasisSet = calc["OptBas"] if jobtype.upper() == 'OPT' else calc["SPBas"]
        TDict   = TheoryDict[LevelTheory]
        BDict   = BasisDict[BasisSet]

        for t in [AllDict,RemDict,TDict,BDict]:
            qcopts['rem'].update(t['rem'])
            if 'basis' in t:
                qcopts['basis'] = t['basis']
            
        qcopts["comments"] = [M.comms[0],LevelTheory+'/'+BasisSet]
        qcopts["molecule"] = ["%i %i" % (calc['Charge'],calc['Mult'])]
        for i in range(M.na()):
            qcopts["molecule"].append(molecule.format_xyz_coord(M.elem[i],M.xyzs[0][i]))

        qcopts["opt"] = ["CONSTRAINT"]
        for i,j in enumerate(val):
            # Unfortunately OpenBabel betrays me here, and it generates dihedrals using the wrong atom indices
            # This is a workaround, read the dihedral back in using the geometry.
            phi = get_dihedral(M.xyzs[0],self.Scans[i][0])
            if n.abs(phi-j) > 1e-2:
                print "Ahooga, OpenBabel created dihedral angle % .4f but really we have % .4f" % (j, phi)
                print "We're just keeping OpenBabel honest.. you may use the OB#|QC# syntax in Control.txt"
                raw_input()
            atoms = " ".join(["%i" % int(k.split('|')[-1]) for k in self.Scans[i][0]])
            qcopts["opt"].append("tors %s %.1f" % (atoms,phi))
        qcopts["opt"].append("ENDCONSTRAINT")

        qcprint = [qc_section(qcopts,"comments"),
                   qc_section(qcopts,"molecule"),
                   qc_rem_section(qcopts)]
        if jobtype == 'opt':
                   qcprint.append(qc_section(qcopts,"opt"))

        if qc_section(qcopts,"basis") != []:
            qcprint.append(qc_section(qcopts,"basis"))

        with open(outfnm,'w') as f: f.write('\n\n'.join(['\n'.join(i) for i in qcprint])+'\n')

    def evaluate(self, calc, val):
        Answer = None
        xyzdir = os.path.join(cwd,'Molecules',self.Dir,*["%+04i" % i for i in val])
        optdir = os.path.join(xyzdir,'%s_%s' % (calc["OptTheory"], calc["OptBas"]))
        spdir = os.path.join(optdir,'%s_%s' % (calc["SPTheory"], calc["SPBas"]))
        mkdir(xyzdir)
        if not os.path.exists(os.path.join(xyzdir,'start.xyz')):
            shutil.copyfile(self.initxyz,os.path.join(cwd,'temp','temp0.xyz'))
            for i,j in enumerate(val):
                atoms = " ".join(["%i" % int(k.split('|')[0]) for k in self.Scans[i][0]])
                o,_ = sp.Popen("%s/obrotate \"%s\" %s/temp/temp%i.xyz %s %i" % (OBPath,self.Smiles,cwd,i,atoms,j),shell=True,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
                with open(os.path.join(cwd,'temp','temp%i.xyz' % (i+1)),'w') as f: f.write(o)
            with open(os.path.join(cwd,'Molecules',self.Dir,'scan_start.xyz'),'a') as f: f.write(o)
            with open(os.path.join(xyzdir,'start.xyz'),'w') as f: f.write(o)
        mkdir(optdir)
        # Check if optimization output already exists
        optout = os.path.join(optdir,'opt.out')
        optstat,optcyc = check_finished(optout)
        optxyz = "This_Is_Not_A_File"
        if optstat == 2:
            print optout, "exists but job is incomplete!"
            touch(os.path.join(optdir,'.FAILED'))
            remove(os.path.join(optdir,'.RUNNING'))
            #raw_input()
        elif optstat == 1: pass # The job appears to be running but the file does not exist (if it's still running, it'll do this)
        elif optstat == 0:
            # If the optimization is done, save the optimized XYZ to a file.
            M = molecule.Molecule(os.path.join(xyzdir,'start.xyz'))
            M.xyzs = [qcparse.get_xyz(optout)[-1]]
            optxyz = os.path.join(optdir,'opt.xyz')
            M.write(fnm=optxyz)
            # Figure out how many 
            touch(os.path.join(optdir,'.COMPLETE'))
            remove(os.path.join(optdir,'.RUNNING'))
        with open('uffe.txt','a') as f: f.write('%i % .3e\n' % ( optcyc, evaluate_uff(os.path.join(xyzdir,'start.xyz'))))
        
        spout = os.path.join(spdir,'sp.out')
        spstat, _ = check_finished(spout)
        if spstat == 2:
            print spout, "exists but job is incomplete!"
            #raw_input()
        elif spstat == 1: pass
        elif spstat == 0:
            Extracted_Data = TheoryDict[calc["SPTheory"]]["get_energy"](spout)
            return Extracted_Data[0]
        #if Extracted_Data[1] == [True]:
        #        return Extracted_Data[0][0]
            
            touch(os.path.join(spdir,'.COMPLETE'))
            remove(os.path.join(spdir,'.RUNNING'))
            #print spdir.strip()
        #print optdir
        #print 
        #if optstat == 0 and spdir == optdir:
        #    print val, qcparse.get_energy_opt(optout)[-1]

        # Create the input file if it doesn't exist
        optin = os.path.join(optdir,'opt.in')
        #preopt = os.path.join(xyzdir,'%s_%s' % ('WB97X-D','6-311+GD'),'opt.xyz')
        
        if not os.path.exists(optin):
            xyzinit = os.path.join(xyzdir,'start.xyz')
            M = molecule.Molecule(xyzinit)
            self.create_qcin(M,calc,val,optin,'opt')
        if os.path.exists(os.path.join(optdir,'.COMPLETE')): pass
        elif os.path.exists(os.path.join(optdir,'.RUNNING')): pass
        elif os.path.exists(os.path.join(optdir,'.FAILED')): pass
        else:
            print "Submitting optimization job in %s" % optdir
            submit_job(optdir,'opt')

        mkdir(spdir)
        # Create the input file if it doesn't exist
        sp_in = os.path.join(spdir,'sp.in')
        if os.path.exists(optxyz):
            if not os.path.exists(sp_in):
                M1 = molecule.Molecule(optxyz)
                self.create_qcin(M1,calc,val,sp_in,'sp')
            if os.path.exists(os.path.join(spdir,'.COMPLETE')): pass
            elif os.path.exists(os.path.join(spdir,'.RUNNING')): pass
            elif os.path.exists(os.path.join(optdir,'.FAILED')): pass
            else:
                print "Submitting single-point job in %s" % spdir
                submit_job(spdir,'sp')

def evaluate_uff(absfnm):
    o,_ = sp.Popen("%s/obenergy -ff UFF %s" % (OBPath,absfnm),shell=True,stdout=sp.PIPE,stderr=sp.PIPE).communicate()
    for line in o.split('\n'):
        if 'TOTAL VAN DER WAALS ENERGY' in line:
            return float(line.split()[-2])
    return "SNOO!"

def check_finished(fnm):
    status =  1
    optcyc = -1
    if os.path.exists(fnm):
        if any(["Thank you very much" in l for l in open(fnm)]): 
            status = 0
        else:
            status = 2
        optcyc = qcparse.get_opt_cycles(fnm)[0]
    return status, optcyc

def submit_job(dnm,fnm):
    ### This block submits the jobs
    os.chdir(dnm)
    task = work_queue.Task('/home/leeping/opt/bin/qchem40 -np 8 %s.in %s.out' % (fnm,fnm))
    task.specify_input_file(os.path.join(dnm,'%s.in' % fnm),'%s.in' % (fnm))
    task.specify_output_file(os.path.join(dnm,'%s.out' % fnm),'%s.out' % (fnm))
    task.specify_algorithm(work_queue.WORK_QUEUE_SCHEDULE_FCFS)
    task.specify_tag(("%s." % fnm) + dnm)
    #print "I'm in submit_job, but skipping"
    touch(os.path.join(dnm,'.RUNNING'))
    wq.submit(task)
    os.chdir(cwd)

def pause(str):
    print str
    raw_input()

def evaluate():
    return 0

def parse(infnm):
    JobList = []
    for line in open(infnm).readlines():
        if re.match('^#',line):
            continue
        elif re.match('^[0-9]+ ',line): # Matches a new simulation line
            S = ScanJob(line)
            JobList.append(S)
        elif re.match('.*//',line): # Matches a calculation line
            S.add_calc(line)
    return JobList

def mkdir(absdir):
    if not os.path.exists(absdir):
        os.makedirs(absdir)

def touch(absfnm):
    if os.path.exists(absfnm):
        os.remove(absfnm)
    f = open(absfnm,'w')
    f.close()

def remove(absfnm):
    if os.path.exists(absfnm):
        os.remove(absfnm)

def wq_wait():
    while not wq.empty():
        print '---'
        task = wq.wait(10)
        if task:
            print 'A job has finished!'
            print 'Job name = ', task.tag, 'command = ', task.command
            print 'output', task.output,
            print 'id', task.id
            print "preferred_host = ", task.preferred_host, 
            print "status = ", task.status, 
            print "return_status = ", task.return_status, 
            print "result = ", task.result, 
            print "host = ", task.host
            #print "submit_time = ", task.submit_time, 
            #print "start_time = ", task.start_time, 
            #print "finish_time = ", task.finish_time
            #print "transfer_start_time = ", task.transfer_start_time, 
            print "computation_time = ", task.computation_time/1000000, 
            print "total_bytes_transferred = ", task.total_bytes_transferred,
            #print "total_transfer_time = ", task.total_transfer_time
            if task.result != 0:
                wq.submit(task)
            else:
                del task
        print "Workers: %i init, %i ready, %i busy, %i total joined, %i total removed" \
            % (wq.stats.workers_init, wq.stats.workers_ready, wq.stats.workers_busy, wq.stats.total_workers_joined, wq.stats.total_workers_removed)
        print "Tasks: %i running, %i waiting, %i total dispatched, %i total complete" \
            % (wq.stats.tasks_running,wq.stats.tasks_waiting,wq.stats.total_tasks_dispatched,wq.stats.total_tasks_complete)
        print "Data: %i / %i kb sent/received" % (wq.stats.total_bytes_sent/1000, wq.stats.total_bytes_received/1024)

def get_dihedral(xyz,atoms):
    r2d = 180./n.pi
    an = [int(k.split('|')[-1])-1 for k in atoms]
    #print an
    #an = [i-1 for i in atoms]
    r1 = xyz[an[1]] - xyz[an[0]]
    r2 = xyz[an[2]] - xyz[an[1]]
    r3 = xyz[an[3]] - xyz[an[2]]
    t1 = n.linalg.norm(r2)*n.dot(r1,n.cross(r2,r3))
    t2 = n.dot(n.cross(r1,r2),n.cross(r2,r3))
    phi = n.arctan2(t1,t2)*r2d
    if phi < -179:
        phi += 360
    return phi

def main():
    JL = parse('Control.txt')
    #JL[0].scan()
    for J in JL:
        J.scan()

    print('THE PORT IS %d' % wq.port)
    wq_wait()


if __name__ == "__main__":
    main()
