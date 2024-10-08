#!/usr/bin/env python

# Imports
#--------
import os
import subprocess
import sys
import time

# Functions
#----------

def makedir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def write_jscript(hpix,marker,dir):
    job_name = "cf_sparse_"+str(hpix)
    job_file = dir+"/jobfiles/"+job_name+".sh"
    # The following code writes lines to the "job_name.sh" file.
    # Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
    # All other lines are read by the shell
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name="+job_name+"\n")                         # job name
        fh.writelines("#SBATCH --output="+dir+"/outfiles/"+job_name+".out\n")      # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+dir+"/outfiles/"+job_name+".err\n")       # error file
        fh.writelines("#SBATCH --partition=unsafe\n")            # queue partition to run the job in
        fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
        fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
        fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem=60000\n")                   # memory in  MB
        fh.writelines("#SBATCH --time=00:60:00\n")               # Maximum job run time
        fh.writelines("module load Anaconda2020.07\n")	         # load anaconda, needed for running python on Hyalite!
        fh.writelines("source activate $HOME/condaenv/\n")       # load conda environment
        fh.writelines("python "+dir+"/files/canfind_v2.py "+hpix+" "+marker)       # write python command
    return job_file,job_name


# Main Code
#----------
if __name__=="__main__":

    hpix = str(sys.argv[1])
    marker = str(sys.argv[2])
    basedir = "/home/x25h971/canfind_dr2"
    makedir(basedir+"/hpix/hgroup_"+str(int(hpix)//1000))
    job_script,job_name = write_jscript(hpix,marker,basedir)
    os.system("sbatch %s" %job_script)
    eflag=0
    while eflag==0:
        print("sleepin for 60")
        time.sleep(60)
        jstat = subprocess.getoutput("sacct -n -X --format state --name "+job_name).split("\n")[-1]
        print("job is "+jstat.strip())
        if jstat!="RUNNING" and jstat!="PENDING" and jstat!="REQUEUED":
            info = subprocess.getoutput("sacct -n -P --delimiter=',' --format jobid,maxrss,maxvmsize,cputimeraw --name "+job_name).split("\n")[-1]
            print("job info = "+info)
            eflag=1
