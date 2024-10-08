#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from astropy import utils, io
from astropy.io import fits
from astropy.table import Table, vstack, join, Column
from dl import queryClient as qc
import healpy as hp
import itertools as it
import math as mat
import math as m # girl what????????????????????????????????????????????
import matplotlib
import numpy as np
from numpy import arange,array,ones,linalg
from sklearn.cluster import DBSCAN
from sklearn import linear_model, datasets
from statistics import mean
import sys

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
def query_fix(query,dtypes):
    '''fix a query in csv format
    Arguments:
       query = the query (str)
       dtypes = data types for each column queried (array-typr)
    '''
    bc=[]
    bd=query.splitlines()
    for j in bd:
        bc.append(j.split(","))
    return(Table(np.array(bc[1:]),names=np.array(bc[0]),dtype=dtypes)) #'f8' makes it a float64

def pm_pass(cluster,out_table):
    ''' Gives object a pass label if it passes the pm test:
    '''
    for kt in range(0,len(cluster)): #for each instance of object
        clonck=out_table['id']==cluster['id'][kt] #using measid, find the same point in out_table
        it=clonck.tolist().index(True) #index of good point
        out_table[it]['pm_test']=1 #give pass label

# Remove a cluster if invalid
def removal(cluster,out_table):
    for kt in range(0,len(cluster)): #for each point in the cluster
        clonck=out_table['measid']==cluster['measid'][kt] #using measid, find the same point in out_table
        it=clonck.tolist().index(True) #index of bad point
        out_table[it]['cluster_label']=-1 #give outlier label

# Give points a track  label 
def labeling(label,out_table,members,pos): #label = some number for the track label, members = track_members function outpit, pos = "p" for pred_pos, "h" for hyp_pos
    for kt in range(0,len(members)): #for each point
        clonck=out_table['measid']==members[kt] #find matching point in out_table
        it=clonck.tolist().index(True) #index of bad point
        if pos=='h': #if hyp_pos was used to predict position,
            out_table[it]['track_h']=label #give appropriate track label 
        if pos=='p':  #if pred_pos was used to predict position,  
            out_table[it]['track_p']=label #give appropriate track label

# Give points a label 
def labelz(members,out_table,column_name,label): #label = some number for the track label, members = track_members function outpit, pos = "p" for pred_pos, "h" for hyp_pos
    for kt in range(0,len(members)): #for each point
        clonck=out_table['measid']==members['measid'][kt] #find matching point in out_table
        it=clonck.tolist().index(True) #index of bad point
        out_table[it][column_name]=label #give appropriate track label 

# RANSAC analysis to remove linearity outliers of cluster in ra-mjd space, 
# also calculates velocity of cluster (at this point, a tracklet) 
def ranslap(cluster):     
    ra=np.array(cluster['ra']) #RA coord.s of cluster members 
    dec=np.array(cluster['dec']) #Dec coords of cluster members
    Xra=np.reshape(ra,(len(ra),1)) #reshape RA coord.s array (X_ra)
    Xdec=np.reshape(dec,(len(dec),1)) #reshape Dec coord.s array  (X_dec)
    y=np.array(cluster['mjd']) #time of measurement of cluster members (Y)

    if len(np.unique(y))<3: #if there are fewer than 3 unique measurement times (fewer than 3 moving object detections)
        inlier_mask=np.zeros(len(ra), dtype=bool) #list of zeros
        outlier_mask=np.ones(len(ra), dtype=bool) #list of ones
        clu=cluster[inlier_mask] #cluster inliers 
        no=cluster[outlier_mask] #cluster outliers
        return(clu,no,0,0) #does not RANSAC but gives all members "outlier" labels
    else: #if there are at least 3 independent moving object detections in the cluster,
        # Robustly fit linear model with RANSAC algorithm-----------------------------------------
        ransac_ra = linear_model.RANSACRegressor(residual_threshold=.0001)
        ransac_ra.fit(Xra, y) #RANSAC on RA and MJD
        ransac_dec = linear_model.RANSACRegressor(residual_threshold=.0001)
        ransac_dec.fit(Xdec, y) #RANSAC on Dec and MJD
        # Predict data of estimated models--------------------------------------------------------
        line_Xra = np.reshape(np.arange(Xra.min(), Xra.max()+.0001,step=((Xra.max()+.0001)-Xra.min())/20),(20,1))
        line_y_ransac_ra = ransac_ra.predict(line_Xra) #line for RANSAC fit, ra
        line_Xdec = np.reshape(np.arange(Xdec.min(), Xdec.max()+.0001,step=((Xdec.max()+.0001)-Xdec.min())/20),(20,1))
        line_y_ransac_dec = ransac_dec.predict(line_Xdec) #line for RANSAC fit, dec
        xsra=np.concatenate(line_Xra) #x values of RA,MJD RANSAC line
        xsdec=np.concatenate(line_Xdec) #x values of Dec,MJD RANSAC line
        ysra=line_y_ransac_ra #y values of RA,MJD RANSAC line
        ysdec=line_y_ransac_dec #y values of Dec,MJD RANSAC line
        mra = (ysra[-1]-ysra[0])/(xsra[-1]-xsra[0]) # 1/slope of RA,MJD RANSAC (1/velocity in RA)
        mdec = (ysdec[-1]-ysdec[0])/(xsdec[-1]-xsdec[0]) # 1/slope of Dec,MJD RANSAC (1/velocity in Dec)
        inlier_mask = ransac_ra.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)
        clu=cluster[inlier_mask] #cluster after outlier removal (cluster inliers)
        no=cluster[outlier_mask] #cluster outliers
        if mra!=0 and mdec!=0: #if the slopes are not zero,
            return (clu,no,1/mra,1/mdec)  #This is the velocity of the object that the tracklet (cluster) represents in ra & dec directions
        else: #if the slopes ARE zero,
            return(clu,no,0,0) #return "0" as velocities

# Validate ``SO tracklets''
def so_check(obj_table,meas_table):  #Get meas table info for each object, determine invalid points in each object and remove, via RANSAC 
    ph=np.array([0])
    ps=np.array([str(0)])
    pf=np.array([float(0)])
    om=0 #counter for tracklet number
    for i in obj_table: #for each object i, except the outliers
        ob_id=i['id']
        if i['pm_test']==1:#if the object passed the pm test, get the object's measurements!
            ms=meas_table['objectid']==ob_id
            if len(meas_table[ms])>2: #if this object more than two mearuements, which should be true from query,
                time=ranslap(meas_table[ms]) #ransac the object measurements!
                if len(time[1])>0: #If there are ransac outliers,
                    removal(time[1],meas_table) #get rid of them
                #Tracklets need to have 3 or more measurements after validation! If so, give them their appropriate velocities.
                cg=meas_table[ms]['cluster_label']!=-1 #$$$$$$$$
                v_ra=i['pmra'] #define tracklet velocity in RA (mas/yr)
                v_dec=i['pmdec'] #define tracklet velocity in dec (mas/yr)
                #print(v_ra,v_dec,type(v_ra),type(v_dec))
                if (len(meas_table[ms][cg])>2) and (v_ra!=float("inf")) and (v_dec!=float("inf")):
                    ##print(np.sqrt(v_ra**2+v_dec**2))
                    ##print(v_ra,v_dec)
                     #Add their cluster labels from db_2
                    labelz(meas_table[ms][cg],meas_table,'cluster_label',om)
                    #meas_table['cluster_label'] #add RA velocity column to out_table
                    labelz(meas_table[ms][cg],meas_table,'v_ra',v_ra) #add Dec velocity column to out_table
                    labelz(meas_table[ms][cg],meas_table,'v_dec',v_dec) #add Dec velocity column to out_table
                    om+=1
                else: removal(meas_table[ms],meas_table)

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # - Setup -
    # get hpix data
    pix = sys.argv[1] #healpix number
    crds=hp.pix2ang(128,int(pix),lonlat=True) # identify healpix coordinates
    RA = crds[0]
    DEC = crds[1]
    # identify NSIDE - default is 128
    nside = 128
    possible_nsides = [128,256] # pix and ring256 columns
    if len(sys.argv)>3:
        nside = int(sys.argv[3])
        if nside not in possible_nsides: nside = 128
    # prepare output directory, filename, and format HP query condition
    outdir="/home/x25h971/canfind_dr2/hpix"
    if nside == 128:
        subdir = int(int(pix)//1000) #assign subdirectory (hpix#/1000)
        outfile="/hgroup_%d/healpix01so_%d.fits" % (int(subdir),int(pix))
        fpix_query = "obj.pix="+pix #format the query condition on object.pix
    elif nside == 256:
        mpix = hp.ang2pix(128,RA,DEC,lonlat=True) #get nside128 pix value
        subdir = int(int(mpix)//1000) #assign subdirectory (hpix#/1000)
        outfile="/hgroup_%d/healpix01so_%d_%d.fits" % (int(subdir),int(mpix),int(pi$
        fpix_query = "obj.ring256="+pix #format the query condition on object.r$

    # - Query -
    qc.set_timeout_request(1800)
    # query the NSC object table (nphot>2) for fpix ring256
    cda = qc.query(sql="".join(["select pmra,pmraerr,pmdec,pmdecerr,id from nsc_dr1.object where nphot>2 and "+fpix_query]),fmt='csv')
    cdat = query_fix(cda,['f8','f8','f8','f8','U19'])
    cdat.add_column(Column(np.zeros(len(cdat),dtype=float)),name="pm_test") #0 for doesn't pass the pm test, 1 for does
    # query the NSC measurements table (nphot>2) for fpix ring256
    cme=qc.query(sql="".join(["select meas.mjd,meas.ra,meas.raerr,meas.dec,meas.decerr,meas.measid,meas.objectid,meas.mag_auto from nsc_dr1.meas as meas JOIN nsc_dr1.object as obj on meas.objectid=obj.id WHERE obj.nphot>2 and "+fpix_query]),fmt='csv')
    cmeas = query_fix(cme,['f8','f8','f8','f8','f8','U18','U19','f8'])
    cmeas.add_column(Column(np.repeat(-99.99,len(cmeas)),name="cluster_label"))
    cmeas.add_column(Column(np.repeat(0,len(cmeas)),name="v_ra")) #add RA velocity column to out_table
    cmeas.add_column(Column(np.repeat(0,len(cmeas)),name="v_dec")) #add Dec velocity column to out_table
    if len(cdat)>0: #if there is actually any data,

#--------------------
# SO check
#--------------------
        for obj in np.unique(cdat['id']): #for every object
            ob=cdat['id']==obj
            #run the pm test!
            pm_ave=np.sqrt((np.mean(np.array(cdat[ob]['pmra']))**2)+(np.mean(np.array((cdat[ob]['pmdec']))**2))) #mas/yr
            if pm_ave>0.0000275*3600*1000*365: #if the object proper motion is greater than the cutoff (0.0000275deg=0.0000275*3600*1000*365mas/yr)
                pm_pass(cdat[ob],cdat) #give object pass labels

#---------------------
# Tracklet Validation
#---------------------
        so_check(cdat,cmeas) #get measurements for objects and RANSAC them!
        if len(np.unique(cmeas['cluster_label']))-1>0: #if there is at least 1 tracklet with 3 mmts,
            cmeas.add_column(Column(np.repeat(pix,len(cmeas))),name='pix') #give them a healpix label
            point_less=cmeas['cluster_label']!=-99.99
            out_less=cmeas[point_less]['cluster_label']!=-1
            cmeas[point_less][out_less].write(outdir+outfile,overwrite=True) #write a fits file with all measurements not associated with SOs, and their info
