#!/usr/bin/env python

# AUTHORS: David Nidever 
#          david.nidever@montana.edu
#          Katie Fasbender
#          katiefasbender@montana.edu

# nsc_meas_wrapper.py will run the NSC Measurements process on all exposures in a given list on tempest.montana.edu


# This script writes a "job_name.sh" file for x exposures, y times every z seconds.
# What that means:
# - In each "job_name.sh" file, "x" exposures will be analyzed
# - "y" number of job files will be submitted in a batch before sleeping for "z" seconds.
# x = 1, y = numjobs, z = 60



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


def write_jscript(job_name,partition,cmd,exp_cmd,dir):
    '''writes a SLURM job script to "job_name.sh"
    Arguments:
    ----------
    job_name (str)
            name of job, job script file
    partition (str)
            node/partition the job will run on
    cmd (str)
            python command to run exposure
    dir (str)
            base directory
    Returns:
    --------
    job_file (str)
            job filename the job script is written to
    '''
    job_file = dir+"outfiles/"+job_name+".sh"
    # The following code writes lines to the "job_name.sh" file.
    # Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
    # All other lines are read by the shell
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        if partition=="priority": fh.writelines("#SBATCH --account=priority-davidnidever\n")       #specify the account to use
        fh.writelines("#SBATCH --job-name="+job_name+"\n")        # job name
        fh.writelines("#SBATCH --output="+dir+"outfiles/"+job_name+".out\n")      # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+dir+"outfiles/"+job_name+".err\n")       # error file
        fh.writelines("#SBATCH --partition="+partition+"\n")     # queue partition to run the job in
        fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
        fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
        fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care
        fh.writelines("#SBATCH --mem=3000\n")                    # memory, set --mem with care!!!!! refer to hyalite quickstart guide
        fh.writelines("#SBATCH --time=48:00:00\n")               # Maximum job run time
        fh.writelines("module load Anaconda3/2021.05\n")         # load anaconda, needed for running python on Hyalite!
        fh.writelines("module load GCC\n")
        fh.writelines("source activate $HOME/condaenv/\n")
        fh.writelines(exp_cmd+"\n")                                   # write python command to download exposure files
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
    parser = ArgumentParser(description='Run NSC Instcal Measurement.')
    parser.add_argument('--version', type=str, nargs=1, help='Version number')
    parser.add_argument('--partition',type=str,nargs=1,help='Delimited list of partitions to divide jobs between')
    parser.add_argument('-r','--redo', action='store_true', help='Redo exposure that were previously processed')
    parser.add_argument('--maxjobs', type=int, nargs=1, default=1, help='Max number of jobs to maintain at any given time, check every 300 sec')
    parser.add_argument('--list',type=str,nargs=1,default=None,help='Input list of exposures to use')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    # Inputs
    version = dln.first_el(args.version)
    redo = args.redo
    print(redo,type(redo))
    ##nmulti = dln.first_el(args.nmulti)
    partitions=args.partition[0].split(',')
    npar=len(partitions)
    maxjobs = int(args.maxjobs[0])//npar
    inputlist = args.list
    if inputlist is not None:
        inputlist = inputlist[0]
    nside = 128
    radeg = 180 / np.pi
    t0 = time.time()

    # Establish necessary directories - figure out for tempest
    basedir = "/home/x25h971/nsc/instcal/"+version+"/"
    mssdir = basedir+"exposures/" #outdir for get_data() func (where exposures will be downloaded)
    localdir = "/home/x25h971/nsc/instcal/"+version+"/"
    outfiledir = basedir+"outfiles/"
    tmpdir = localdir+"tmp/"

    makedir(mssdir)
    makedir(outfiledir)
    makedir(tmpdir)
    subdirs = ['logs','c4d','k4m','ksb']
    for sub in subdirs:
        makedir(basedir+sub)

    # Log File 
    #---------
    # Create Log file name;
    # format is nsc_combine_main.DATETIME.log
    ltime = time.localtime()
    # time.struct_time(tm_year=2019, tm_mon=7, tm_mday=22, tm_hour=0, tm_min=30, tm_sec=20, tm_wday=0, tm_yday=203, tm_isdst=1)
    smonth = str(ltime[1])
    if ltime[1]<10: smonth = '0'+smonth
    sday = str(ltime[2])
    if ltime[2]<10: sday = '0'+sday
    syear = str(ltime[0])[2:]
    shour = str(ltime[3])
    if ltime[3]<10: shour='0'+shour
    sminute = str(ltime[4])
    if ltime[4]<10: sminute='0'+sminute
    ssecond = str(int(ltime[5]))
    if ltime[5]<10: ssecond='0'+ssecond
    logtime = smonth+sday+syear+shour+sminute+ssecond
    logfile = basedir+'logs/nsc_instcal_measure_main.'+logtime+'.log'

    # Set up logging to screen and logfile
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()
    fileHandler = logging.FileHandler(logfile)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)
    rootLogger.setLevel(logging.NOTSET)

    rootLogger.info("Creating NOAO InstCal Measurement catalogs")
    rootLogger.info("host = "+host)
    rootLogger.info("version = "+version)
    ##rootLogger.info("nmulti = "+str(nmulti))
    ##rootLogger.info("hosts = "+','.join(np.atleast_1d(hosts)))
    rootLogger.info("partitions = "+str(partitions))
    rootLogger.info("redo = "+str(redo))

    # Loading the exposure list(s)
    #-----------------------------
    if inputlist is None:
        rootLogger.info('Using input lists:')
        rootLogger.info('  '+basedir+'/lists/decam_instcal_list.fits.gz')
        rootLogger.info('  '+basedir+'/lists/mosaic3_instcal_list.fits.gz')
        rootLogger.info('  '+basedir+'/lists/bok90prime_instcal_list.fits.gz')
        list1 = fits.getdata(basedir+'/lists/decam_instcal_list.fits.gz',1)
        list2 = fits.getdata(basedir+'/lists/mosaic3_instcal_list.fits.gz',1)
        list3 = fits.getdata(basedir+'/lists/bok90prime_instcal_list.fits.gz',1)
        lstr = dln.concatenate([list1,list2,list3])    # concatenate
        del(list1,list2,list3)
    else:
        rootLogger.info('Using input list: '+inputlist)
        lstr = fits.getdata(inputlist,1) # r for random (about to do)
        # Check that it has all the columns that we need
        needcols = ['INSTRUMENT','FLUXFILE','MASKFILE','WTFILE','DATE_OBS']
        for n in needcols:
            if n not in lstr.dtype.names:
                raise ValueError('Column '+n+' not in '+inputlist)
    nlstr = dln.size(lstr)
    rootLogger.info(str(nlstr)+' InstCal images')

    # Putting them in RANDOM but REPEATABLE order
    rootLogger.info('RANDOMIZING WITH SEED=1')
    rnd = np.arange(nlstr)
    np.random.seed(1)
    np.random.shuffle(rnd)
    lstr = lstr[rnd]

    gdexp = np.arange(nlstr)
    ngdexp = nlstr


    # Check the exposures
    #--------------------
    rootLogger.info('Checking on the exposures')
    dtype_expstr = np.dtype([('instrument',str,100),('rawname',str,100),('fluxfile',str,100),('wtfile',str,100),
                             ('maskfile',str,100),('outfile',str,150),('done',bool),('torun',bool),('cmd',str,1000),('cmddir',str,1000),('cputime',str,100),
                             ('exp_cmd',str,1000),('submitted',bool),('jobname',str,100),('jobid',str,100),('partition',str,100),('jobstatus',str,100),
                             ('maxrss',str,100),('maxvmsize',str,100)])
    expstr = np.zeros(ngdexp,dtype=dtype_expstr)  # string array for exposure info
    expstr['jobid']="-99.99"
    for i in range(ngdexp):
        #if i % 500 == 0: rootLogger.info(i)

        # Format exposure info for string array
        instrument = lstr['INSTRUMENT'][gdexp[i]].strip()
        if type(instrument) is bytes: instrument=instrument.decode()
        rawname = lstr['RAWNAME'][gdexp[i]].strip()
        if type(rawname) is bytes: rawname=rawname.decode()
        fluxfile = lstr['FLUXFILE'][gdexp[i]].strip()
        if type(fluxfile) is bytes: fluxfile=fluxfile.decode()
        wtfile = lstr['WTFILE'][gdexp[i]].strip()
        if type(wtfile) is bytes: wtfile=wtfile.decode()
        maskfile = lstr['MASKFILE'][gdexp[i]].strip()
        if type(maskfile) is bytes: maskfile=maskfile.decode()
        fdir,base = os.path.split(fluxfile)

        # Change the root directory name to reflect host repo structure
        # format on tempest will be basedir+/exposures/
        ##lo = fluxfile.find('/mss1/')
        ##fluxfile = mssdir+fluxfile[(lo+6):]
        fluxfile = fluxfile.split('/')[-1]
        ##lo = wtfile.find('/mss1/')
        ##wtfile = mssdir+wtfile[(lo+6):]
        wtfile = wtfile.split('/')[-1]
        ##lo = maskfile.find('/mss1/')
        ##maskfile = mssdir+maskfile[(lo+6):]
        maskfile = maskfile.split('/')[-1]

        expstr['instrument'][i] = instrument
        expstr['rawname'][i] = rawname
        expstr['fluxfile'][i] = fluxfile
        expstr['wtfile'][i] = wtfile
        expstr['maskfile'][i] = maskfile

        # Check if the output already exists.
        dateobs = lstr['DATE_OBS'][gdexp[i]]
        if type(dateobs) is np.bytes_: dateobs=dateobs.decode()
        night = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
        baseroot = base[0:base.find('.fits.fz')]
        outfile = basedir+instrument+'/'+night+'/'+baseroot+'/'+baseroot+'_1.fits'
        expstr['outfile'][i] = outfile
        #print("outfile = ",outfile)

        # Does the output file exist?
        expstr['done'][i] = False
        if os.path.exists(outfile): expstr['done'][i] = True
        if expstr['fluxfile'][i]==expstr['wtfile'][i]: expstr['done'][i] = True

        # If no outfile exists or yes redo:
        if (expstr['done'][i]==False) or (redo==True):
            #if file_test(file_dirname(outfile),/directory) eq 0 then file_mkdir,file_dirname(outfile)  ; make directory
            expstr['exp_cmd'][i] = 'python '+basedir+'get_exposures.py '+rawname+' '+fluxfile+' '+wtfile+' '+maskfile+' '+mssdir
            expstr['cmd'][i] = 'python '+basedir+'nsc_instcal_meas.py '+fluxfile+' '+wtfile+' '+maskfile+' '+version
            #expstr['exp_cmd'][i] = 'echo $PATH'
            #expstr['exp_cmd'][i] = 'python /home/x25h971/nsc/instcal/v4/shortrun.py'
            #print(expstr['exp_cmd'][i])
            #print(expstr['cmd'][i])
            #print(rawname,fluxfile,wtfile,maskfile)
            expstr['cmddir'][i] = tmpdir
            expstr['torun'][i] = True
        # If outfile exists and no redo:
        elif (expstr['done'][i]==True) and (redo==False):
            expstr['torun'][i] = False
            #rootLogger.info(outfile+' EXISTS and --redo NOT set')

    # Parcel out jobs
    #----------------
    # Define exposures to run & total #jobs/partition
    torun,nalltorun = dln.where(expstr['torun'] == True)    # Total number of jobs to run (# exposures)
    ntorun = len(torun)
    rootLogger.info(str(ntorun)+" exposures")
    if ntorun == 0:
        rootLogger.info('No exposures to process.')
        sys.exit()
    njpar = 0 #number of jobs per partition (divide evenly)
    if ntorun!=0: njpar = ntorun//(maxjobs*npar)
    rootLogger.info(str(ntorun)+' exposures to process on '+str(maxjobs*npar)+' "partitions", and '+str(npar)+' real partitions.')
    #rootLogger.info('Maintaining '+str(njobs)+' job(s) submitted per partition.')
    sleep_time=10     #"z" # of seconds to sleep between checking on job batches

    # Split exposures evenly among defined partitions
    par_arr = np.array([i*(maxjobs*npar) for i in range(0,ntorun//(maxjobs*npar))]) # array of every 'npar'th index
    partitions = np.reshape([[i+"_"+str(parchan) for i in partitions] for parchan in range(0,maxjobs)],(maxjobs*npar))
    print("partitions = ",partitions)
    for part,i in zip(partitions,range(0,maxjobs*npar)):
        expstr['partition'][torun[par_arr+i]] = part

    #print(expstr)
    runfile = basedir+'/nsc_instcal_measure_main.'+logtime+'_run.fits'
    Table(expstr).write(runfile)
    # Start submitting jobs
    #----------------------
    jb = 0
    while (jb < ntorun):
        for part in partitions:
            rootLogger.info("Checking status of last job submitted to "+part+" partition")
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
                rootLogger.info("Job id="+last_jid+" is still "+lj_status+", sleepin for 30")
                time.sleep(sleep_time)
            # ---If last job is completed, failed, cancelled, killed, or none: submit a new job!
            else:
                rootLogger.info("--Submitting new job to "+part+" partition--")
                # if last job was completed, get some info about it 
                if lj_status=="COMPLETED":
                    ljinfo = subprocess.getoutput("sacct -n -P --delimiter=',' --format cputimeraw,maxrss,maxvmsize --jobs "+last_jid)
                    ljinfo = ljinfo.split("\n")[-1].split(",")
                    expstr['cputime'][torun[lsub]] = ljinfo[0]
                    expstr['maxrss'][torun[lsub]] = ljinfo[1]
                    expstr['maxvmsize'][torun[lsub]] = ljinfo[2]
                    # remove downloaded files!
                    if os.path.exists(mssdir+expstr['fluxfile'][torun[lsub]].strip()) and os.path.exists(mssdir+expstr['wtfile'][torun[lsub]].strip()) and os.path.exists(mssdir+expstr['maskfile'][torun[lsub]].strip()):
                        os.remove(mssdir+expstr['fluxfile'][torun[lsub]].strip())
                        os.remove(mssdir+expstr['wtfile'][torun[lsub]].strip())
                        os.remove(mssdir+expstr['maskfile'][torun[lsub]].strip())

                # get index & info of next job to submit
                next_sub = list(partition_ind & unsubmitted_ind)
                print("length of next_sub array = ",len(next_sub))
                if len(next_sub)==0: jbsub = ntorun-1
                else: jbsub = np.sort(next_sub)[0]

                # check to see if the exposures have been downloaded yet...
                #if not, carry on to the next partion
                # if they have, create and submit the job!
                if os.path.exists(mssdir+expstr['fluxfile'][torun[jbsub]].strip()) and os.path.exists(mssdir+expstr['wtfile'][torun[jbsub]].strip()) and os.path.exists(mssdir+expstr['maskfile'][torun[jbsub]].strip()): 
                    cmd = expstr['cmd'][torun[jbsub]]
                    exp_cmd = expstr['exp_cmd'][torun[jbsub]]
                    cmddir = expstr['cmddir'][torun[jbsub]]
                    partition = expstr['partition'][torun[jbsub]].split("_")[0]

                    # --Write job script to file--
                    job_name = 'nsc_meas_'+str(logtime)+'_'+str(jb)
                    job_file=write_jscript(job_name,partition,cmd,exp_cmd,basedir)

                    # --Submit job to slurm queue--
                    os.system("sbatch %s" %job_file)
                    expstr['submitted'][torun[jbsub]] = True
                    rootLogger.info("Job "+job_name+"  submitted to "+part+" partition; sleeping for 30 sec")
                    time.sleep(sleep_time) #let the job get submitted 
                    # get jobid of submitted job, update array of exposure info
                    jid = subprocess.getoutput("sacct -n -X --format jobid --name "+job_name)
                    jid = jid.split("\n")[-1]
                    expstr['jobname'][torun[jbsub]] = job_name
                    expstr['jobid'][torun[jbsub]] = jid
                    jb+=1

                # save job structure, sleep before checking/submitting again 
                runfile = basedir+'/nsc_instcal_measure_main.'+logtime+'_run.fits'
                if os.path.exists(runfile): os.remove(runfile)
            Table(expstr).write(runfile,overwrite=True)
        #time.sleep(sleep_time)  # SLEEP before checking/submitting next jobs 


