#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# final_concat.py concatenates the hgroup128 tracklet fits files from CANFind
# to create one final tracklet catalog, cfdr<dr#>_tracklet_cat.fits.gz

#--------------------
# Imports
#--------------------
from astropy.table import Table,vstack
import numpy as np
import os

#--------------------
# Main Code
#--------------------

if __name__=="__main__":

    basedir = "/home/x25h971/canfind_dr2/files/"
    concatdir = "/home/x25h971/canfind_dr2/concats/"
    t_mmts = Table()
    t_tlets = Table()
    for root,dir,files in os.walk(concatdir):
        for name in files:
            fname = os.path.join(root,name)
            if name.split("_")[-1].split(".")[0]=="mmts":
                print("adding mmts from ",fname)
                tab = Table.read(fname)
                t_mmts = vstack([t_mmts,tab])
            if name.split("_")[-1].split(".")[0]=="tracklets":
                print("adding tracklets from ",fname)
                tab = Table.read(fname)
                t_tlets = vstack([t_tlets,tab])
    t_mmts.write(concatdir+"cfdr2_tracklet_mmts.fits.gz",overwrite=True)
    print("tracklet measurements ",concatdir+"cfdr2_tracklet_mmts.fits.gz"," written")
    t_tlets.write(concatdir+"cfdr2_tracklet_cat.fits.gz",overwrite=True)
    print("tracklet ",concatdir+"cfdr2_tracklet_cat.fits.gz"," written")
