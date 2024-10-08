#!/usr/bin/env

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# This script will create slurm jobs to create tracklet pair lists.
# Each job will produce one list of tracklet pairs from the NOIRLab Source Catalog (NSC),
# for one HEALPix (NSIDE=32), using the script 'make_tpair_lists.py'

# Input = full list of HEALPix (NSIDE=32) for the NSC, fits file/astropy table

#-------------
# Imports
#-------------
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord
from astropy.table import Column,Row,Table,vstack
from astropy.time import Time
import astropy.units as u
from dlnpyutils import utils as dln
import healpy as hp
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import subprocess
import sys
import time


#-------------
# Functions
#-------------

def makedir(dir):
    '''makes a directory with name "dir" if it does not exist
    Arguments:
    ----------
    dir (str)
            directory name
    Returns:
    --------
        None; directory "dir" is created if not already there
    '''
    if not os.path.exists(dir):
        os.mkdir(dir)

def write_jscript(job_name,cmd,dir):
    '''writes a SLURM job script to "job_name.sh"
       Lines starting with #SBATCH are read by Slurm.
       Lines starting with ## are comments.
       All other lines are read by the shell
    Arguments:
    ----------
    job_name (str)
            name of job, job script file
    cmd (str)
            python command
    dir (str)
            where to write the job script & outfiles
    Returns:
    --------
    job_file (str)
            job filename the job script is written to
    '''
    job_file = dir+"/outfiles/"+job_name+".sh"
    # The following code writes lines to the "job_name.sh" file.
    # Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
    # All other lines are read by the shell
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        #fh.writelines("#SBATCH --account=priority-davidnidever\n")       #specify the account to use
        fh.writelines("#SBATCH --job-name="+job_name+"\n")                         # job name
        fh.writelines("#SBATCH --output="+dir+"outfiles/"+job_name+".out\n")      # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+dir+"outfiles/"+job_name+".err\n")	   # error file
        fh.writelines("#SBATCH --partition=unsafe\n")            # queue partition to run the job in
        fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
        fh.writelines("#SBATCH --cpus-per-task=1\n")             # cpus per task
        fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
        fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem=3000\n")                   # memory in  MB
        fh.writelines("#SBATCH --time=5:00:00\n")               # Maximum job run time
        fh.writelines("module load Anaconda3\n")           # load anaconda, needed for running python on Hyalite!
        fh.writelines("source activate $HOME/condaenv/\n")	 # load conda environment
        fh.writelines(cmd.strip()+'\n')                      # write python CF command
        return job_file



#-------------
# Main Code
#-------------
if __name__=="__main__":

    #Inputs
    # Initiate input arguments
    parser = ArgumentParser(description='Create lists of tracklet pairs for orbit calculation')
    parser.add_argument('--dr', type=int, nargs=1, help='NSC Data Release')
    parser.add_argument('--partition',type=str,nargs=1,help='Delimited list of partitions to divide jobs between')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposure that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=1, help='Max number of jobs to maintain at any given time')
    parser.add_argument('--hp32',type=str,nargs=1,default=None,help='Input list of HEALPix (NSIDE32) to use')
    args = parser.parse_args()

    # Start time
    t0 = time.time()

    # Inputs
    redo = args.redo
    dr = args.dr[0]
    print(redo,type(redo))
    partitions=args.partition[0].split(',')
    npar=len(partitions)
    maxjobs = int(args.maxjobs[0])//npar #max number of jobs to maintain on a partition
    hp32file = args.hp32
    if hp32file is not None:
        hp32file = hp32file[0]
    hp32list = Table.read(hp32file)
    hp32list['done'] = Column(np.repeat(False,len(hp32list)))
    hp32list['torun'] = Column(np.repeat(False,len(hp32list)))
    hp32list['submitted'] = Column(np.repeat(False,len(hp32list)))
    hp32list['cmd'] = Column(dtype="U500",length=len(hp32list))
    hp32list['outfile'] = Column(dtype="U500",length=len(hp32list))
    hp32list['jobname'] = Column(dtype="U200",length=len(hp32list))
    hp32list['jobid'] = Column(np.repeat("-99.99",len(hp32list)))
    hp32list['partition'] = Column(dtype="U50",length=len(hp32list))
    hp32list['jobstatus'] = Column(dtype="U50",length=len(hp32list))


    basedir = "/home/x25h971/"
    hgroup32dir = basedir+"orbits_dr"+str(dr)+"/tpair_lists/hgroup32_" #+str(pix32//1000)+"/" is added later in code
    concatdir = basedir+"canfind_dr2/concats/"
    fbase = concatdir+"cf_dr2_hgroup_" # the base name for tracklet_concat files

    # for each pix32, check completeness & get info
    print("Creating tracklet pair lists for "+str(len(hp32list))+" HEALPix (NSIDE32)")
    for p32 in range(len(hp32list['pix32'])):
        pix32 = hp32list['pix32'][p32]
        outfile = hgroup32dir+str(pix32//1000)+"/cfdr2_"+str(pix32)+"_tpairs.fits"
        if os.path.exists(outfile): hp32list['done'][p32] = True
        hp32list['outfile'][p32] = outfile
        if hp32list['done'][p32]==False or redo==True: # to run
            hp32list['torun'][p32] = True
            hp32list['cmd'][p32] = "python "+basedir+"canfind_dr2/files/make_tpair_list.py "+str(pix32)+" "+str(dr)
        elif hp32list['done'][p32]==True and redo==False: hp32list['torun'][p32] = False


    # Parcel out jobs
    #----------------
    # Define exposures to run & total #jobs/partition
    torun,nalltorun = dln.where(hp32list['torun'] == True)    # Total number of jobs to run (# exposures)
    ntorun = len(torun)
    print(str(ntorun)+" PIX32")
    if ntorun == 0:
        print('No PIX32 left')
        sys.exit()
    njpar = 0 #number of jobs per partition (divide evenly)
    if ntorun!=0: njpar = ntorun//(maxjobs*npar)
    print(str(ntorun)+' tracklet lists to make on '+str(maxjobs*npar)+' "portitions", and '+str(npar)+' real partitions.')
    #rootLogger.info('Maintaining '+str(njobs)+' job(s) submitted per partition.')
    sleep_time=15     #"z" # of seconds to sleep between checking on job batches

    # Split exposures evenly among defined partitions
    par_arr = np.array([i*(maxjobs*npar) for i in range(0,ntorun//(maxjobs*npar))]) # array of every 'npar'th index
    partitions = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,maxjobs)],(maxjobs*npar))
    print("partitions = ",partitions)
    for part,i in zip(partitions,range(0,maxjobs*npar)):
        hp32list['partition'][torun[par_arr+i]] = part

    runfile = basedir+"orbits_dr"+str(dr)+"/tpair_lists/maketpair."+str(t0)+"_run.fits"
    hp32list.write(runfile)
    # Start submitting jobs
    #----------------------
    jb = 0
    while (jb < ntorun):
        for part in partitions:
            print("Checking status of last job submitted to "+part+" partition")
            partition_ind = set(np.where(hp32list['partition'][torun]==part)[0])
            submitted_ind = set(np.where(hp32list['submitted'][torun]==1)[0])
            unsubmitted_ind = set(np.where(hp32list['submitted'][torun]==0)[0])

            # get index & status of last job submitted
            last_sub = list(partition_ind & submitted_ind)
            if len(last_sub)==0:
                lsub=np.sort(list(partition_ind))[0]
            else:
                lsub=np.sort(last_sub)[-1]
            last_jid = hp32list[torun[lsub]]['jobid']
            last_jname = hp32list[torun[lsub]]['jobname']
            if last_jid != "-99.99": lj_status = (subprocess.getoutput("sacct -n -X --format state --jobs="+last_jid).split("\n")[-1]).strip()
            else: lj_status = "NONE" #no jobs have been submitted
            hp32list[torun[lsub]]['jobstatus'] = lj_status
            print("lj_status = ",lj_status,", jobname = ",last_jname,last_jid)
            # ---If last job is still running: wait!
            if (lj_status=="RUNNING" or lj_status=="PENDING" or lj_status=="REQUEUED"):
                print("Job id="+last_jid+" is still "+lj_status+", sleepin for a few")
                time.sleep(sleep_time)
            # ---If last job is completed, failed, cancelled, killed, or none: submit a new job!
            else:
                print("--Submitting new job to "+part+" partition--")
                # if last job was completed, get some info about it 
                if lj_status=="COMPLETED":
                    #ljinfo = subprocess.getoutput("sacct -n -P --delimiter=',' --format cputimeraw,maxrss,maxvmsize --jobs "+last_jid)
                    #ljinfo = ljinfo.split("\n")[-1].split(",")
                    hp32list['done'][torun[lsub]] = True

                # get index & info of next job to submit
                next_sub = list(partition_ind & unsubmitted_ind)
                print("length of next_sub array = ",len(next_sub))
                if len(next_sub)==0: jbsub = ntorun-1
                else: jbsub = np.sort(next_sub)[0]

                # create and submit the job!
                pix32 = hp32list['pix32'][torun[jbsub]]
                cmd = hp32list['cmd'][torun[jbsub]]
                jdir = basedir+"orbits_dr"+str(dr)+"/tpair_lists/"
                partition = hp32list['partition'][torun[jbsub]].split("_")[0]

                # --Write job script to file--
                job_name = "tpairs_"+str(pix32)
                job_file = write_jscript(job_name,cmd,jdir)

                # --Submit job to slurm queue--
                os.system("sbatch %s" %job_file)
                hp32list['submitted'][torun[jbsub]] = True
                print("Job "+job_name+"  submitted to "+part+" partition; sleeping for a few")
                time.sleep(sleep_time) #let the job get submitted 
                # get jobid of submitted job, update array of exposure info
                jid = subprocess.getoutput("sacct -n -X --format jobid --name "+job_name)
                jid = jid.split("\n")[-1].strip()
                hp32list['jobname'][torun[jbsub]] = job_name
                hp32list['jobid'][torun[jbsub]] = jid
                jb+=1
                hp32list.write(runfile,overwrite=True)






