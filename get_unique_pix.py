#!/usr/bin/env python

from astropy.table import Table
import healpy as hp
import numpy as np
import os
import subprocess

if __name__=="__main__":

    # get list of all CF output files (contain pix128 values)
    g=subprocess.getoutput("ls hgroup_* ")
    g=g.split("\n")
    g=np.array(g)
    bools=[g[i][0:6]!="hgroup" for i in range(len(g))]
    g=g[bools]
    g=g[g!=""]

    hplist = Table(names=("pix128","pix64","pix32"),dtype=["int","int","int"])
    #hp128list=[int(g[i].split(".fits")[0].split("_")[-1]) for i in range(len(g))]

    for file in g:
        file_split = file.split(".fits")[0].split("_")
        if len(file_split)==3: pix128 = int(file_split[-1])
        elif len(file_split)==4: pix128 = int(file_split[-2])
        coords = hp.pix2ang(128,pix128,lonlat=True)
        pix64 = hp.ang2pix(64,coords[0],coords[1],lonlat=True)
        pix32 = hp.ang2pix(32,coords[0],coords[1],lonlat=True)
        hplist.add_row([pix128,pix64,pix32])
    print(len(hplist)," unique pix128 values,")
    print(len(np.unique(hplist['pix32']))," unique pix32 values")
    hplist.write("cfdr2_unique_pix.fits",overwrite=True)


