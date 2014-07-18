"""
This module provides API for simple cluster job submission.
The only user currently is genome_validator.py.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import sys
import os
import subprocess
import shutil

class Cluster:
    """Provides methods for creating and submitting cluster jobs"""
    def __init__(self, head_node):
        self.head_node = head_node
        
    def create_job(self, job_file, cmds, name=None, queue='all.q', mem='1G', working_dir=None, depend=[], ncpus=1, excl=False, email=None):
        """Creates job file"""
        if ',' in mem:
            mem, ncpus = mem.split(',')
        
        out = open(job_file, 'w')
        out.write("#! /bin/bash\n")
        out.write("#$ -S /bin/bash\n")
        out.write("#$ -R y\n")
        out.write("#$ -q %s\n" % (queue))
        out.write("#$ -l mem_free=%s -l mem_token=%s -l h_vmem=%s\n" % (mem, mem, mem))
        
        if name:
            out.write("#$ -N %s\n" % (name))
        if working_dir:
            out.write("#$ -wd %s\n" % (working_dir))
        if ncpus is not None and ncpus > 1:
            out.write("#$ -pe ncpus %s\n" % (ncpus))
        if excl:
            out.write("#$ -l excl=true\n")
        if depend:
            out.write("#$ -hold_jid %s\n" % (','.join([str(jid) for jid in depend])))
        if email:
            out.write("#$ -M %s\n#$ -m as\n" % (email))
            
        out.write('\n')
        for cmd in cmds:
            out.write(cmd + '\n')
        out.close()
        
        return True
         
    def submit_job(self, job_file):
        """Submits job given job file"""       
        sys.stderr.write('submitting %s to %s\n' % (job_file, self.head_node))
                
        proc = subprocess.Popen(
            ['ssh', self.head_node, 'qsub', job_file],
            shell=False, cwd=os.getcwd(), stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            )
        proc.wait()
        
        submit_stdout = proc.stdout.readlines()
        submit_stderr = proc.stderr.readlines()
        
        if submit_stdout:
            for mess in submit_stdout:
                print mess

        if submit_stderr:
            for mess in submit_stderr:
                print mess
        
        job_id = None
        if "Your job" in submit_stdout[0]:
            job_id = int(submit_stdout[0].split()[2])

        return job_id
        