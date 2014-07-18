#!/usr/bin/env python

# written by:    Ka Ming Nip
# updated on:    July 11, 2014

import string
import re
import os
import sys
import subprocess
import stat
from utilities import cfg

class MyJob:    
    def __init__(self, content=None, id=None, name=None, predecessors=None, config_key=None, threads=1, first_task_id=None, last_task_id=None, queue=None, outputfiles=None):
        self.name = name
        self.id = id
        self.content = content
        self.predecessors = predecessors
        self.config_key = config_key
        self.threads = threads
        self.first_task_id = first_task_id
        self.last_task_id = last_task_id
        self.queue = queue
        self.outputfiles = outputfiles


def getJobScriptConfigs(config_file):
    configs = {}
    fh = open(config_file, 'r')
    for line in fh:
        line_stripped = line.strip()
        if len(line_stripped) > 0:
            if ':' in line_stripped:
                items = line_stripped.split(':', 1)
                configs[items[0].strip()] = items[1].strip()
    fh.close()
    
    #some sanity checks
    for key in ["local", "cluster_basic", "cluster_parallel", "cluster_basic_array", "cluster_parallel_array", "predecessors_list_delimiter", "submit_cluster_job_return_string", "submit_cluster_array_job_return_string", "submit_cluster_job_command", "run_local_job_command", "cluster_max_resources", "cluster_max_resources_array"]:
        if not configs.has_key(key) or not configs[key]:
            sys.stderr.write("ERROR: \"%s\" is not defined in %s" % (key, config_file))
            sys.exit(1)
    
    if not "${JOBID}" in configs["submit_cluster_job_return_string"]:
        sys.stderr.write("ERROR: Cannot find \"${JOBID}\" in the \"submit_cluster_job_return_string\" defined in " + config_file) 
        sys.exit(1)
    
    if not "${JOBID}" in configs["submit_cluster_array_job_return_string"]:
        sys.stderr.write("ERROR: Cannot find \"${JOBID}\" in the \"submit_cluster_array_job_return_string\" defined in " + config_file) 
        sys.exit(1)
    
    return configs


def makeJobScript(template, content, jobname=None, workingdir=None, threads=1, predecessors=None, memory=None, setup=None, first_task_id=None, last_task_id=None, tmpmem=None, queue=None, logdir=None, error_exit_status=100):
    script = ""
    if not os.path.isfile(template):
        sys.stderr.write("ERROR: Cannot find template file, " + template)
        sys.exit(1)
    fh = open(template, 'r')
    use_wd = False
    for line in fh:
        #Replace known variable names in template with the appropriate values        
        if "${JOB_NAME}" in line:
            if jobname:
                line = line.replace("${JOB_NAME}", str(jobname))
            else:                
                continue
        if "${WORKING_DIR}" in line:
            if workingdir:
                line = line.replace("${WORKING_DIR}", str(workingdir))
                use_wd = True
            else:
                continue
        if "${THREADS}" in line:
            if threads and threads > 1:
                line = line.replace("${THREADS}", str(threads))
            else:
                continue
        if "${PREDECESSORS}" in line:
            if predecessors:
                line = line.replace("${PREDECESSORS}", str(predecessors))
            else:
                continue
        if "${MEM}" in line:
            if memory:
                line = line.replace("${MEM}", str(memory))
            else:
                continue
        if "${SETUP_PATHS}" in line:
            if setup:
                line = line.replace("${SETUP_PATHS}", str(setup))
            else:
                continue
        if "${FIRST_TASK_ID}" in line:
            if first_task_id:
                line = line.replace("${FIRST_TASK_ID}", str(first_task_id))
            else:
                continue
        if "${LAST_TASK_ID}" in line:
            if last_task_id:
                line = line.replace("${LAST_TASK_ID}", str(last_task_id))
            else:
                continue
        if "${QUEUE}" in line:
            if queue:
                line = line.replace("${QUEUE}", str(queue))
            else:
                continue
        if "${CONTENT}" in line:            
            if content:
                content_str = '{ ' + ' && '.join(content) + '; } || { echo ERROR: command exit status $? ; exit %d; }' % error_exit_status
                line = line.replace("${CONTENT}", content_str)
            else:
                continue
        if "${TMPMEM}" in line:
            if tmpmem:
                line = line.replace("${TMPMEM}", str(tmpmem))
            else:
                continue
        if "${LOG_DIR}" in line:
            if logdir:
                if workingdir and use_wd:
                    if logdir.rstrip(os.sep) == workingdir.rstrip(os.sep):
                        continue
                    if logdir.startswith(workingdir):
                        logdir = logdir.replace(workingdir, '', 1)
                        logdir = logdir.lstrip(os.sep)
                line = line.replace("${LOG_DIR}", str(logdir))
            else:
                continue
        script += line
    fh.close()
    return script


def run_locally(jobs, job_script_config_file, job_script_template_dir, prefix=None, logdir=None, debug=False, workingdir=None):
    job_script_configs = getJobScriptConfigs(job_script_config_file)
    
    print "Combining jobs for local run ..."
    template = os.path.join(job_script_template_dir, job_script_configs["local"])
    
    content = []
    for j in jobs:        
        if j.first_task_id and j.last_task_id:
            tmp_cmd_list = []            
            tmp_cmd_list.extend(j.content)
            tmp_cmd_list[0] = "for TA_JOBID in {%s..%s};do " % (j.first_task_id, j.last_task_id) + tmp_cmd_list[0]
            tmp_cmd_list[-1] = tmp_cmd_list[-1] + "; done"
            content.extend(tmp_cmd_list)
        else:
            content.extend(j.content)
            
    script = makeJobScript(template=template, setup=None, content=content, workingdir=workingdir)
    
    name = "local"
    if prefix != None:
        name = prefix + ".local"
    
    scriptpath = '%s/%s.sh' % (logdir, name)
    print "script: %s" % scriptpath
    
    # write script to file
    fh = open(scriptpath, 'w')
    fh.write(script)
    fh.close()
    
    permission = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR # owner read/write/execute
    permission = permission | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP # group read/write/execute
    os.chmod(scriptpath, permission)
    
    cmd = [job_script_configs["run_local_job_command"], scriptpath]
    
    # run script locally
    if debug:
        print " ".join(cmd)
    else:
        print "Run script locally ..."
        exit_status = subprocess.call(cmd)
        if exit_status != 0:
            sys.stderr.write("ERROR: Execution of script ended with a non-zero exit-status.")
            sys.exit(1)
        print "Finished running script."


def submit_jobs(jobs, settings, job_script_config_file, job_script_template_dir, head_node=None, logdir=None, debug=False, queue=None, workingdir=None, max_resources=False):
    job_script_configs = getJobScriptConfigs(job_script_config_file)

    predecessors_list_delimiter = str(job_script_configs["predecessors_list_delimiter"])
    
    submit_cluster_return_string = job_script_configs["submit_cluster_job_return_string"]
    cluster_basic_template_file = os.path.join(job_script_template_dir, job_script_configs["cluster_basic"])
    cluster_parallel_template_file = os.path.join(job_script_template_dir, job_script_configs["cluster_parallel"])

    submit_cluster_return_string_array = job_script_configs["submit_cluster_array_job_return_string"]
    cluster_basic_template_file_array = os.path.join(job_script_template_dir, job_script_configs["cluster_basic_array"])
    cluster_parallel_template_file_array = os.path.join(job_script_template_dir, job_script_configs["cluster_parallel_array"])
    
    cluster_max_template_file = os.path.join(job_script_template_dir, job_script_configs["cluster_max_resources"])
    cluster_max_template_file_array = os.path.join(job_script_template_dir, job_script_configs["cluster_max_resources_array"])

    submit_cluster_job_command = job_script_configs["submit_cluster_job_command"]
    error_exit_status = int(job_script_configs["job_error_status"])
    
    prog = re.compile("\s*" + submit_cluster_return_string.replace("${JOBID}", "(?P<jid>\S+)") + "\s*")
    prog_array = re.compile("\s*" + submit_cluster_return_string_array.replace("${JOBID}", "(?P<jid>\w+)") + "\s*")

    for job in jobs:
        predecessors_string = None
        if job.predecessors and len(job.predecessors) > 0:
            predecessor_ids = []
            for j in job.predecessors:
                if j and j.id:
                    predecessor_ids.append(str(j.id))
            if len(predecessor_ids) > 0:
                predecessors_string = predecessors_list_delimiter.join(predecessor_ids)
                            
        memory = None
        tmpmem = None
        threads = job.threads
        if job.config_key:
            memory = cfg.get_value(settings, 'memory', job.config_key)
            tmpmem = cfg.get_value(settings, 'tmpmem', job.config_key)
            if memory:
                items = memory.split(',')
                num_items = len(items)
                # memory
                if num_items > 0:
                    memory_str_stripped = items[0].strip()
                    if len(memory_str_stripped) > 0:
                        memory = memory_str_stripped
                # threads
                if num_items > 1:
                    threads_str_stripped = items[1].strip()
                    if len(threads_str_stripped) > 0:
                        threads = int(threads_str_stripped)
        
        template = cluster_basic_template_file
        matcher = prog
        
        if max_resources or (memory is not None and memory.lower() == 'max'):
            template = cluster_max_template_file
        elif threads > 1:
            template = cluster_parallel_template_file
        #endif
        
        if job.first_task_id and job.last_task_id:
            matcher = prog_array
            template = cluster_basic_template_file_array
            if max_resources or (memory is not None and memory.lower() == 'max'):
                template = cluster_max_template_file_array
            elif threads > 1:
                template = cluster_parallel_template_file_array
            #endif
        #endif
        
        if job.queue:
            queue = job.queue
        #endif
                
        #template, content, jobname=None, workingdir=None, threads=1, predecessors=None, memory=None, setup=None
        script = makeJobScript(template=template, content=job.content, jobname=job.name, logdir=logdir, threads=threads, predecessors=predecessors_string, first_task_id=job.first_task_id, last_task_id=job.last_task_id, memory=memory, setup=None, queue=queue, tmpmem=tmpmem, workingdir=workingdir, error_exit_status=error_exit_status)
        scriptpath = logdir + os.sep + job.name + ".sh"
        print "job script: " + scriptpath
        fh = open(scriptpath, 'w')
        fh.write(script)
        fh.close()
        
        permission = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR # owner read/write/execute
        permission = permission | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP # group read/write/execute
        os.chmod(scriptpath, permission)
        
        cmd = []
        submit_options = cfg.get_value(settings, 'commands', submit_cluster_job_command).split()
        if head_node:
            #ssh to head node then submit the job with qsub            
            cmd = ["ssh", head_node, submit_cluster_job_command]
        else:
            #submit the job with qsub
            cmd = [submit_cluster_job_command]
        
        if len(submit_options) > 0:
            cmd.extend(submit_options)
            
        cmd.append(scriptpath)
                
        print " ".join(cmd)
        
        if not debug:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout, stderr = p.communicate()
            
            jobid = None
            
            lines = stdout.split("\n")            
            for line in lines:
                print line
                m = matcher.match(line)
                if m:
                    jobid = m.group('jid')
            
            if not jobid:
                sys.stderr.write("ERROR: Cannot identify job id.\n")
                sys.exit(1) 
            
            job.id = jobid

