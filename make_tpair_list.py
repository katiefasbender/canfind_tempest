#!/usr/bin/env

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu


# This script will create a list of tracklet pairs for one HEALPix (NSIDE=32)
# from the NOIRLab Source Catalog (NSC).
# Input = one HEALPix (NSIDE=32) from the NSC (pix32)
# Output = full list of tracklet pairs for the pix32

#-------------
# Imports
#-------------
from astropy.coordinates import SkyCoord
from astropy.table import Column,Row,Table,vstack
from astropy.time import Time
import astropy.units as u
import healpy as hp
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
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

def read_cfmpc80(filename):
    '''Read a text file (filename.txt) written in "CANFind" MPC 80-column format, 
       which means the first 13 columns are CF-style measurement ids.  
    Arguments:
    ----------
    filename (str)
        Name of file to be read, <filename>.txt
    Returns:
    --------
    mmt_cat (astropy table)
        Astropy table with 1 row for each mmt in the file, and the following columns:
        - mmt_id (measurement id) XX.YYY.ZZZZZZ
        - tracklet_id (tracklet id)  YYY.ZZZZZZ where
              X = unique mmt # within tracklet "Y"
              Y = unique tracklet # within HP "Z"
              Z = unique HEALPIX # (NSIDE=64, nested ordering)
        - line (corresponding line without measurement/tracklet ids, for Find_Orb purposes)
        - mjd, ra, dec, mag_auto, filter
    '''
    # Make sure the file exists
    fil=filename
    if os.path.exists(fil) is False:
        print(fil+" NOT found")
        return None
    lines = readlines(fil)
    nmmts = len(lines)-2
    if nmmts == 0:
        print("No measurements in "+file)
        return None
    mmt_cat=Table(names=("mmt_id","tracklet_id","line","mjd","ra","dec","mag_augo","filter"),
                  dtype=["U13","U10","U92","float64","float64","float64","float64","U1"]) # output table (cat)
    for ln in lines[2:]: 
        mmt_id = ln[0:13] # get the measurement(observation) id 
        tracklet_id = mmt_id[3:13]
        #tracklet_id_ln = "   "+ln[3:92]
        temp_id_ln = ln[13:] #"            "
        # get string info
        yyyy=ln[15:19]
        mm=ln[20:22]
        dd=ln[23:31]
        ra_hh=ln[32:34]
        ra_mm=ln[35:37]
        ra_ss=ln[38:44]
        dec_sign=ln[44]
        dec_dd=ln[45:47]
        dec_mm=ln[48:50]
        dec_ss=ln[51:56]
        mag=ln[65:70]
        filt=ln[70]
        ra_err=ln[81:86]
        dec_err=ln[87:92]
        # transform strings to values
        ddd = float(dd) - int(float(dd))
        h = int(ddd*24)
        m = int(((ddd*24)-h)*60)
        s = ((((ddd*24)-h)*60)-m)*60
        mjd = Time({'year': int(yyyy), 'month': int(mm), 'day': int(float(dd)),'hour': h, 'minute': m, 'second': s},scale='utc')
        mjd = mjd.mjd
        c = SkyCoord(ra_hh+'h'+ra_mm+'m'+ra_ss+'s', dec_sign+dec_dd+'d'+dec_mm+'m'+dec_ss+'s', frame='icrs')
        ra= c.ra.value
        dec=c.dec.value
        mag_val=float(mag)
        # add row to table
        row=[mmt_id,tracklet_id,temp_id_ln,mjd,ra,dec,mag_val,filt]
        mmt_cat.add_row(row)
    return(mmt_cat)


#-------------
# Main Code
#-------------
if __name__=="__main__":

    #Inputs
    dr = sys.argv[2]           # NSC data release
    pix32 = int(sys.argv[1])   # HP32 (nested ordering)
    print("making tracklet pair list for pix32 = ",pix32)

    basedir = "/home/x25h971/"
    hgroup32dir = basedir+"orbits_dr"+str(dr)+"/tpair_lists/hgroup32_" #+str(pix32//1000)+"/" is added later in code
    concatdir = basedir+"canfind_dr2/concats/"
    fbase = concatdir+"cf_dr2_hgroup_" # the base name for tracklet_concat files

    # for the pix32, make a list of tracklet pairs to fit to orbits
    tpid = 0 #a unique number for each testable tracklet pair in the pix32 (fo_id, aka pair id)
    tracklet_pairs = Table(names=("tracklet_0","tracklet_1","hgroup_filename0","hgroup_filename1","pix32","fo_id","ntracklets","mjd","dmjd"),
                            dtype=["U10","U10","U150","U150","float64","U100","int32","float64","float64"]) # a table to hold the tracklet pairs list


    # --get the corresponding 16 NSIDE128 hpix--
    coords32 = hp.pix2ang(32,int(pix32),lonlat=True)  # get the coordinates for the nside=32 pix
    nbs128=hp.get_all_neighbours(128,coords32[0],coords32[1],lonlat=True)     # get the 8 nearest neighbors to the cooordinates for nside=128    
    coords128 = hp.pix2ang(128,nbs128,lonlat=True)    # get the center coordinates for the 8 nside=128 neighbors
    pix64 = np.unique(hp.ang2pix(64,coords128[0],coords128[1],lonlat=True))   # find the 4 unique corresponding nside=64 healpix
    coords64=hp.pix2ang(64,pix64,lonlat=True)         # get the center coordinates for the 4 nside=64  neighbors
    nbs256=hp.get_all_neighbours(256,coords64[0],coords64[1],lonlat=True)     # get the 8*4=32 nearest neighbors to the cooordinates for nside=256
    coords256=hp.pix2ang(256,nbs256,lonlat=True)      # get the center coordinates for the 32 nside=256  neighbors
    pix128s=np.unique(hp.ang2pix(128,coords256[0],coords256[1],lonlat=True))  # find the 16 unique corresponding nside=128 helapix 

    print("pix128 values = ",pix128s)
    # --get the list of hgroup128s for this pix32--
    hgroups128=np.unique(pix128s//1000)

    # --get all mmts & tracklets for this pix32--
    mmts_32 = Table()        # a table for the pix32 tracklet measurements
    tlets_32 = Table()       # a table for the pix32 tracklets
    for hg in hgroups128:    #for each pix128 in this pix32, get the mmts and tracklets 
        if os.path.exists(fbase+str(hg)+"_tracklets.fits") and os.path.exists(fbase+str(hg)+"_mmts.fits"):
            print("reading tracklet & mmt files for hgroup128 = ",hg)
            tlets_hgroup128=Table.read(fbase+str(hg)+"_tracklets.fits") #read in the tracklet & mmt concat files 
            mmts_hgroup128=Table.read(fbase+str(hg)+"_mmts.fits")
            # get the mmts/tracklets with the correct hpix
            #print(tlets_hgroup128['pix128'])
            # mmts
            mbools = [(mmts_hgroup128['pix128']==pix128s[i]) for i in range(len(pix128s))]
            mpx = np.repeat(False,len(mmts_hgroup128))
            for mb in range(len(mbools)):
                mpx = np.logical_or(mpx,mbools[mb])
            # tlets
            tbools = [(tlets_hgroup128['pix128']==pix128s[i]) for i in range(len(pix128s))]
            tpx = np.repeat(False,len(tlets_hgroup128))
            for tb in range(len(tbools)):
                tpx = np.logical_or(tpx,tbools[tb])
            mmts_32=vstack([mmts_32,mmts_hgroup128[mpx]])
            tlets_32=vstack([tlets_32,tlets_hgroup128[tpx]])

    # --make list of all tracklet pair combinations--
    if len(tlets_32)>1:
        print(len(tlets_32)," tracklets")
        tpairs=list(combinations(tlets_32['tracklet_id'],2))

        # --prune the list--
        for pr in tpairs:
            # get info for each tracklet
            t0,t1=pr #name each tracklet, t0 and t1
            mms0=np.where(mmts_32['tracklet_id']==t0)[0]
            mms1=np.where(mmts_32['tracklet_id']==t1)[0]
            mms=np.concatenate([mms0,mms1])
            t_mmts = mmts_32[mms]
            t_mmts.sort('mjd')
            tinfo0=np.where(tlets_32['tracklet_id']==t0)[0]
            tinfo1=np.where(tlets_32['tracklet_id']==t1)[0]
            tinfo=np.concatenate([tinfo0,tinfo1])
            tlets=tlets_32[tinfo]
            tlets.sort('mjd')
            # --prune 1----------------------------------------------------------------------------------------------
            # --each mmt of both tracklets must have a unique date
            cond1 = (len(np.unique(t_mmts['mjd']))==len(t_mmts))
            # --prune 2------------------------------------------------------------------------------------------
            # --the first mmt date of the tracklet that was detected second must be >= the last mmt date of the first-detected tracklet 
            # --(mjd_2>=mjd_1+dmjd+1)
            cond2 = ((tlets[0]['mjd']+tlets[0]['dmjd'])<tlets[1]['mjd'])
            if cond1 and cond2: #ok! this pair can be tested for an orbit
                hp_subdirs = [int(t[-6:])//1000 for t in pr]
                fname0 = fbase+str(hp_subdirs[0])+".txt"
                fname1 = fbase+str(hp_subdirs[1])+".txt"
                t_mjd0 = t_mmts['mjd'][0]
                t_dmjd = t_mmts['mjd'][-1]-t_mjd0
                tracklet_pairs.add_row([tlets['tracklet_id'][0],tlets['tracklet_id'][1],fname0,fname1,pix32,str(tpid),2,t_mjd0,t_dmjd])
                tpid+=1
        makedir(hgroup32dir+str(pix32//1000))
        tlist_name = hgroup32dir+str(pix32//1000)+"/cfdr2_"+str(pix32)+"_tpairs.fits"
        tracklet_pairs.write(tlist_name,overwrite=True)
        print("tracklet pair list "+tlist_name+" written for pix32 = ",str(pix32))
