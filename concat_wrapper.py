#!/usr/bin/env python

# AUTHOR:  Katie Fasbender
#          katiefasbender@montana.edu

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from argparse import ArgumentParser
from astropy.table import Table,Column
from astropy.io import fits
from dlnpyutils import utils as dln, coords
import logging
import numpy as np
import os
import socket
import subprocess
import sys
import time

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

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


def write_jscript(job_name,partition,cmd,outdir):
    '''writes a SLURM job script to "job_name.sh"
    Arguments:
    ----------
    job_name (str)
            name of job, job script file
    partition (str)
            node/partition the job will run on
    cmd (str)
            python command to run exposure
    outdir (str)
            base directory
    Returns:
    --------
    job_file (str)
            job filename the job script is written to
    '''
    job_file = outdir+job_name+".sh"
    # The following code writes lines to the "job_name.sh" file.
    # Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
    # All other lines are read by the shell
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        if partition=="priority": fh.writelines("#SBATCH --account=priority-davidnidever\n")       #specify the account to use
        fh.writelines("#SBATCH --job-name="+job_name+"\n")        # job name
        fh.writelines("#SBATCH --output="+outdir+job_name+".out\n")      # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+outdir+job_name+".err\n")       # error file
        fh.writelines("#SBATCH --partition="+partition+"\n")     # queue partition to run the job in
        fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
        fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
        fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem=6000\n")                    # memory, set --mem with care!!!!! refer to hyalite quickstart guide
        fh.writelines("#SBATCH --time=6:00:00\n")               # Maximum job run time
        fh.writelines("module load Anaconda3\n")         # load anaconda, needed for running python on Hyalite!
        fh.writelines("source activate $HOME/condaenv/\n")
        fh.writelines(cmd+"\n")                                       # write python command to analyze exposure
        fh.writelines("conda deactivate")
    return job_file




#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Concatenate CANFind tracklet files.')
    parser.add_argument('--partition',type=str,nargs=1,help='Delimited list of partitions to divide jobs between')
    parser.add_argument('-r','--redo', action='store_true', help='Redo concats that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=1, help='Max number of jobs to maintain at any given time, check every 300 sec')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list hgroup128 (pix128//1000)')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()

    # Inputs
    partitions=args.partition[0].split(',')
    redo=args.redo
    print("redo = ",redo)
    npar=len(partitions)
    maxjobs = int(args.maxjobs[0])//npar
    inputlist = args.list
    if inputlist is not None:
        inputlist = inputlist[0]

    # Establish necessary directories - figure out for tempest
    basedir = "/home/x25h971/canfind_dr2/"
    localdir = basedir+"files/"
    outdir = basedir+"concats/"
    outfiledir = outdir+"outfiles/"

    makedir(outfiledir)

    # Load the list of hgroups
    hgroup_list = Table.read(inputlist)
    nhgroup = len(hgroup_list)

    # Prepare hgroup job info
    #--------------------
    dtype_expstr = np.dtype([('hgroup',str,100),('done',bool),('torun',bool),('cmd',str,1000),('cmddir',str,1000),
                             ('submitted',bool),('jobname',str,100),('jobid',str,100),('partition',str,100),
                             ('jobstatus',str,100),('cputime',str,100),('maxrss',str,100),('maxvmsize',str,100),
                             ('outfile_mpc80',str,100),('outfile_t',str,100),('outfile_m',str,100)])
    expstr = np.zeros(nhgroup,dtype=dtype_expstr)  # string array for exposure info
    expstr['jobid']="-99.99"
    for i in range(nhgroup):


        expstr['hgroup'][i] = hgroup_list['hgroup'][i]
        expstr['outfile_mpc80'][i] = outdir+"cf_dr2_hgroup_"+str(hgroup_list['hgroup'][i])+".txt" #the MPC 80-col format output file
        expstr['outfile_m'][i] = outdir+"cf_dr2_hgroup_"+str(hgroup_list['hgroup'][i])+"_mmts.fits" #the MPC 80-col format output file
        expstr['outfile_t'][i] = outdir+"cf_dr2_hgroup_"+str(hgroup_list['hgroup'][i])+"_tracklets.fits" #the MPC 80-col format output file
        outfiles = [expstr['outfile_mpc80'][i],expstr['outfile_m'][i],expstr['outfile_t'][i]]
        # Check if the output already exists.
        expstr['done'][i] = False
        outbools = [os.path.exists(ofile) for ofile in outfiles]
        if False not in outbools: expstr['done'][i] = True

        # If no outfile exists or yes redo:
        if (expstr['done'][i]==False) or (redo==True):
            expstr['cmd'][i] = 'python '+localdir+'tracklet_concat.py '+str(hgroup_list['hgroup'][i])+" 2"
            expstr['cmddir'][i] = localdir
            expstr['torun'][i] = True
        # If outfile exists and no redo:
        elif (expstr['done'][i]==True) and (redo==False):
            expstr['torun'][i] = False

    # Parcel out jobs
    #----------------
    # Define exposures to run & total #jobs/partition
    torun,nalltorun = dln.where(expstr['torun'] == True)    # Total number of jobs to run (# exposures)
    ntorun = len(torun)
    print("running tracklet_concat.py on "+str(ntorun)+" hgroups")
    if ntorun == 0:
        print("no more hgroups to run!")
        sys.exit()
    njpar = 0 #number of jobs per partition (divide evenly)
    if ntorun!=0: njpar = ntorun//(maxjobs*npar)
    print(str(ntorun)+' hgroups to process on '+str(maxjobs*npar)+' "partitions", and '+str(npar)+' real partitions.')
    sleep_time=10     #"z" # of seconds to sleep between checking on job batches

    # Split exposures evenly among defined partitions
    par_arr = np.array([i*(maxjobs*npar) for i in range(0,ntorun//(maxjobs*npar))]) # array of every 'npar'th index
    partitions = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,maxjobs)],(maxjobs*npar))
    print("partitions = ",partitions)
    for part,i in zip(partitions,range(0,maxjobs*npar)):
        expstr['partition'][torun[par_arr+i]] = part

    #print(expstr)
    runfile = localdir+'cfdr2_trackletconcat'+str(t0)+'_run.fits'
    Table(expstr).write(runfile)
    # Start submitting jobs
    #----------------------
    jb = 0
    while (jb < ntorun):
        for part in partitions:
            print("Checking status of last job submitted to "+part+" partition")
            partition_ind = set(np.where(expstr['partition'][torun]==part)[0])
            submitted_ind = set(np.where(expstr['submitted'][torun]==1)[0])
            unsubmitted_ind = set(np.where(expstr['submitted'][torun]==0)[0])

            # get index & status of last job submitted
            last_sub = list(partition_ind & submitted_ind)
            if len(last_sub)==0:
                lsub=np.sort(list(partition_ind))[0]
            else:
                lsub=np.sort(last_sub)[-1]
            last_jid = expstr[torun[lsub]]['jobid']
            last_jname = expstr[torun[lsub]]['jobname']
            if last_jid != "-99.99": lj_status = (subprocess.getoutput("sacct -n -X --format state --jobs="+last_jid).split("\n")[-1]).strip()
            else: lj_status = "NONE" #no jobs have been submitted
            expstr[torun[lsub]]['jobstatus'] = lj_status
            print("lj_status = ",lj_status,", jobname = ",last_jname,last_jid)
            # ---If last job is still running: wait!
            if (lj_status=="RUNNING" or lj_status=="PENDING" or lj_status=="REQUEUED"):
                print("Job id="+last_jid+" is still "+lj_status+", sleepin for a sec")
                time.sleep(sleep_time)
            # ---If last job is completed, failed, cancelled, killed, or none: submit a new job!
            else:
                print("--Submitting new job to "+part+" partition--")
                # if last job was completed, get some info about it 
                if lj_status=="COMPLETED":
                    ljinfo = subprocess.getoutput("sacct -n -P --delimiter=',' --format cputimeraw,maxrss,maxvmsize --jobs "+last_jid)
                    ljinfo = ljinfo.split("\n")[-1].split(",")
                    expstr['cputime'][torun[lsub]] = ljinfo[0]
                    expstr['maxrss'][torun[lsub]] = ljinfo[1]
                    expstr['maxvmsize'][torun[lsub]] = ljinfo[2]

                # get index & info of next job to submit
                next_sub = list(partition_ind & unsubmitted_ind)
                print("length of next_sub array = ",len(next_sub))
                if len(next_sub)==0: jbsub = ntorun-1
                else: jbsub = np.sort(next_sub)[0]

                # check to see if the exposures have been downloaded yet...
                #if not, carry on to the next partion
                # if they have, create and submit the job!
                cmd = expstr['cmd'][torun[jbsub]]
                partition = expstr['partition'][torun[jbsub]].split("_")[0]

                # --Write job script to file--
                job_name = 'cfdr2_concat_'+str(t0)+'_'+str(jb)
                job_file=write_jscript(job_name,partition,cmd,outfiledir)

                # --Submit job to slurm queue--
                os.system("sbatch %s" %job_file)
                expstr['submitted'][torun[jbsub]] = True
                print("Job "+job_name+"  submitted to "+part+" partition; sleeping for a sec")
                time.sleep(sleep_time) #let the job get submitted 
                # get jobid of submitted job, update array of exposure info
                jid = subprocess.getoutput("sacct -n -X --format jobid --name "+job_name)
                jid = jid.split("\n")[-1]
                expstr['jobname'][torun[jbsub]] = job_name
                expstr['jobid'][torun[jbsub]] = jid
                jb+=1

                # save job structure, sleep before checking/submitting again 
            if os.path.exists(runfile): os.remove(runfile)
            Table(expstr).write(runfile,overwrite=True)
        #time.sleep(sleep_time)  # SLEEP before checking/submitting next jobs 


