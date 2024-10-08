#!/usr/bin/env python

# Author: Katie Fasbender
#         katiefasbender@montana.edu

# cf_makecat.py reads in the CANFind tracklet measurements file,
# cfdr2_tracklet_mmts.fits.gz, and writes them to individual
# HEALPix (NSIDE=RING32) files, creating a catalog sorted by ring32.

#-----------
# Imports
#-----------
from astropy.table import Table,Column,vstack
import healpy as hp
import numpy as np
import sys
import os

#-----------
# Functions
#-----------
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


#-----------
# Main Code
#-----------

if __name__=="__main__":

    basedir = "/home/x25h971/"
    cfdir = basedir+"canfind_dr2/files/"
    outdir = basedir+"catalogs/canfind/dr2/ring32/"

    cf_filename = cfdir+"cfdr2_tracklet_mmts.fits.gz"
    cf_mmts = Table.read(cf_filename)
    cf_mmts['mpc_match'] = Column(length=len(cf_mmts),dtype="U100")
    cf_mmts['ring32'] = Column(hp.ang2pix(32,cf_mmts['ra'],cf_mmts['dec'],lonlat=True,nest=False))
    r32s = np.unique(cf_mmts['ring32'])
    print("writing ",len(r32s)," HP RING32 files for ",len(cf_mmts)," CANFind tracklet mmts")
    for r32 in r32s:
        subdir = outdir+str(r32//1000)+"/"
        makedir(subdir)
        r32_cf_filename = subdir+str(r32)+".fits.gz"
        r32_mmts = cf_mmts[cf_mmts['ring32']==r32]
        r32_mmts.write(r32_cf_filename,overwrite=True)
        print("mmts written to ",r32_cf_filename)
