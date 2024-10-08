#!usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# trans_from_hyalite.py is a script that will, for every hgroup_subdir,
# move all accessible CANFind_dr1 healpix files on Hyalite to Tempest.

# COMMAND: $python trans_to_tempest.py </path/to/HP_list>.fits

#---------------------------
# Imports
#---------------------------
from astropy.table import Table
import numpy as np
import os
import sys
import time

#---------------------------
# Functions
#---------------------------

def write_jscript(job_name,cmd,partition,dir):
    '''Writes a job script file "job_name.sh" with given command to given partition
    Arguments:
    ----------
    job_name (str)
        self-explanatory; if not, contact me immediately
    cmd (str)
        command to be written in job script file
    partition (str)
        the hyalite partition to which we plan to submit this job
    dir (str)
        the directory to put the outfiles & job scripts in
    Returns:
    --------
    job_file (str)
        name of job script file just written
    '''
    job_file=dir+job_name+".sh"
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name="+job_name+"\n")              # job name, pix is the #_of_healpix
        fh.writelines("#SBATCH --output="+dir+job_name+".out\n")            # standard output file (%j = jobid, %x should be job name?)
        fh.writelines("#SBATCH --error="+dir+job_name+".err\n")             # standard error file
        fh.writelines("#SBATCH --partition="+partition+"\n")	 # queue partition to run the job in
        fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel? no...
        fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
        fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem=1000\n")                    # MB of Memory allocated; set --mem with care
        fh.writelines("#SBATCH --time=24:00:00\n")               # Maximum job run time
        fh.writelines(cmd)
    return(job_file)


#---------------------------
# Main Code
#---------------------------

if __name__=="__main__":

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------    # --Set up definitions and shit--
    hp_listfile = sys.argv[1]
    hp_list = Table.read(hp_listfile)
    subdir_list = np.unique(hp_list['pix']//1000)
    tempest_basedir = "/home/x25h971/canfind_dr1/hpix/" #where the hgroup_subdirs will go
    tempest_outdir = "/home/x25h971/canfind_dr1/outfiles/" #where job outfiles & scriptfiles will be written
    hyalite_basedir = "/mnt/lustrefs/scratch/katie.fasbender/canfind_hp/canfind_hpix/" #where the hgroup_subdirs are


    # --Loop through each unique subdirectory--
    # to submit a job to transfer each healpix file available
    for sbdr in subdir_list:
        partition = "unsafe"
        hpdir = "hgroup_"+str(sbdr)+"/"
        cmd = "scp -r x25h971@hyalite.rci.montana.edu:"+hyalite_basedir+hpdir+" "+tempest_basedir
        job_name = "hgroup_"+str(sbdr)+"_transfer"
        job_file = write_jscript(job_name,cmd,partition,tempest_outdir)
        os.system("sbatch "+job_file)
        print("submitting job "+job_name+", sleeping")
        time.sleep(30)



