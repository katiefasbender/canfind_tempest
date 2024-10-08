#!/usr/bin/env python
from astropy.table import Table,Column,Row,vstack
import healpy as hp
import numpy as np
import sys

if __name__=="__main__":

    basedir = "/home/x25h971/canfind_dr2/files/"
    hp_file = sys.argv[1]
    hp_base = hp_file.split(".")[0]
    hp_outfile = hp_base+"_split.fits"
    hp_list = Table.read(hp_file)
    nhp128 = len(hp_list) # total number of HP (NSIDE=128)
    hp_list['pix256'] = Column(np.zeros(shape=nhp128))

    print("Splitting up dense HP...")
    # split each dense HP (NMEAS>1e6) into its 4 ring256 counterparts
    job_list = Table()
    for hp128 in range(0,nhp128):
        row128 = hp_list[hp128]
        if hp_list[hp128]['NMEAS']>1000000:
            hp_list[hp128]['pix256']==1
            coords128 = hp.pix2ang(128,hp_list[hp128]['PIX'],lonlat=True)
            neighbors512 = hp.get_all_neighbours(512,coords128[0],coords128[1],lonlat=True) #get the 8 nearest neighbors to the cooordinate$
            coords512 = hp.pix2ang(512,neighbors512,lonlat=True) #get the center coordinates for the 8 nside=512-neighbors
            pix256 = np.unique(hp.ang2pix(256,coords512[0],coords512[1],lonlat=True)) #find the 4 unique corresponding nside=256-hpix
            #uh now add new rows
            for hp256 in pix256:
                row128['pix256'] = int(hp256)
                job_list = vstack([job_list,row128])
            if np.mod(hp128,1000)==0: print("hpix number = ",hp128)
    job_list = vstack([job_list,hp_list[hp_list['pix256']==0]])
    print("New HP file = "+str(len(job_list))+" rows")
    job_list.write(basedir+hp_outfile,overwrite=True)
    print("File written to "+basedir+hp_outfile)
