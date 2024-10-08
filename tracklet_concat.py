#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# tracklet_concat.py is a script that reads fits files in a given subdir (hgroup_<subdir#>) with CANFind output,
# where each fits file contains the tracklet information from one HEALPIX (NSIDE=128).  



# grabs the tracklets, and records them in the Minor Planet Center's 80-column format,
# explained in no detail below (find more at https://www.minorplanetcenter.net/iau/info/OpticalObs.html),
# then writes the new tracklet table to a text file.  This is the file type that digest2 and MPC require.
# It also queries any necessary information? We'll see about that. The tracklet measurements are also
# written to a FITS file, to later be concatenated into the tracklet catalog by cf_tracklet_cat_creator.py

# command format:  $python /path/to/tracklet_concat.py <subdir#> <dr#>

# output files:     1) cf_dr<#>_<subir#>_mmts.fits            one line for each tracklet mmt (cat_mmts)
# (written to       2) cf_dr<#>_<subdir#>_tracklets.fits      one line for each tracklet (cat_tracklets) (still must do)
# basedir/concats/) 3) cf_dr<#>_<subdir#>_mmts_mpc80.txt      same as 1, in MPC 80-col format (mpc_mmts)


# From the MPC (that website listed above):
#---------------------------
#   Columns     Format   Use
#---------------------------
#   For Minor Planets:   (assume this for all tracklets, initially )
#    1 -  5       A5     Packed minor planet number
#    6 - 12       A7     Packed provisional designation, or a temporary designation
#    13           A1     Discovery asterisk
#   For Comets:
#    1 -  4       I4     Periodic comet number
#    5            A1     Letter indicating type of orbit
#    6 - 12       A7     Provisional or temporary designation
#    13           X      Not used, must be blank
#   For Natural Satellites:
#    1            A1     Planet identifier [Only if numbered]
#    2 -  4       I3     Satellite number  [Only if numbered]
#    5            A1     "S"
#    6 - 12       A7     Provisional or temporary designation
#    13           X      Not used, must be blank
#---------------------------
#   Columns     Format   Use
#---------------------------
#   For all: (columns 14-80)
#    14            A1    Note 1
#    15            A1    Note 2
#    16 - 32             Date of observation
#    33 - 44             Observed RA (J2000.0) in deg
#    45 - 56             Observed Dec (J2000.0) in deg
#    57 - 65       9X    Must be blank
#    66 - 71    F5.2,A1  Observed magnitude and band
#    72 - 77       X     Must be blank
#    78 - 80       A3    Observatory code


# Tholen-style uncertainty addition for Find_Orb (MPC will NOT accept these columns!)
#---------------------------
#   Columns     Format   Use
#---------------------------
#    81                  Blank
#    82 - 86             Uncertainty in RA in arcsec
#    87                  Blank
#    88 - 92             Uncertainty in Dec in arcsec



# Example Row:

# 123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|12
# <-ObjDesig->*nnYYYY MM DD.DDDDD HH MM SS.SSSsdd mm ss.ss<blanks >MM.MMBz<ref>COD a.aaa a.aaa

# J009S         C1992 02 24.41844 10 41 36.62 +09 39 12.9          18.6 Vi40908691 0.124 0.026


#---------------------------------------------------------------------------------
# Imports
#---------------------------------------------------------------------------------

import astropy.units as u
from astropy.table import Table,Column,join,vstack
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.time import Time
from dl import queryClient as qc
import healpy as hp
import matplotlib
import numpy as np
from scipy.optimize import curve_fit, least_squares
import sys
import time
import os

#---------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------
def makedir(dir):
    '''Makes a directory with name "dir" if it does not exist.
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

def query_fix(quer,dtps):
    '''Turns a csv-format query into an astropy table
    Arguments:
    ----------
    quer (str)
            the query text (for NSC Measurements table)
    dtps (list)
            desired datatypes for the columns being queried 
    Returns:
    --------
    an astropy table with the queried columns 
    '''
    bc=[]
    bdd=quer.splitlines()
    for ag in bdd:
        bc.append(ag.split(","))
    return(Table(np.array(bc[1:]),names=np.array(bc[0]),dtype=dtps))

def poly(x,coef,*args):
    y = np.array(x).copy()*0.0
    if len(args)>0:
        coef = np.hstack((coef,np.array(args)))
    n = len(coef)
    for ii in range(n):
        y += coef[ii]*x**(n-1-ii)
    return y


def cfit(x,y,nord,sigma=None):
    '''Estimates proper motion in x and y
    Arguments:
    ----------
    Returns:
    --------
    '''
    initpar = np.zeros(nord+1,float)
    bounds=(-np.inf,np.inf)
    try: cfitcoef, cfitcov = curve_fit(poly, x, y, p0=initpar, sigma=sigma, bounds=bounds)
    except: cfitcoef, cfitcov = [np.array([-9999,-9999]),np.array([-9999,-9999])]
    if sigma is not None:
        cfitperr = np.sqrt(np.diag(cfitcov))
        return cfitcoef, cfitperr, cfitcov
    else:return cfitcoef, np.array([-9999,-9999]), cfitcov

def tracklet_calc(tracklet,pix128,ring256,t_num):
    '''Calculate tracklet info
    Arguments:
    ----------
    Returns:
    --------
    '''
# Dates/times
    first_obs_mjd = tracklet[0]['mjd']                          # MJD of 1st tracklet observation [day]
    dt_mjd = tracklet[-1]['mjd'] - tracklet[0]['mjd']
    mean_mjd = np.mean(tracklet['mjd'])
    tracklet_id_val = str(t_num).zfill(3)+"."+str(pix128).zfill(6)
    nmeas = len(tracklet)
# Position in equatorial & ecliptic coords
    mean_ra = np.mean(tracklet['ra'])                           # equatorial tracklet coordinates
    mean_dec = np.mean(tracklet['dec'])
    ring256 = hp.ang2pix(256,mean_ra,mean_dec,lonlat=True,nest=False)
    pix64 = hp.ang2pix(64,mean_ra,mean_dec,lonlat=True)
    pix32 = hp.ang2pix(32,mean_ra,mean_dec,lonlat=True)
    eq_coords=SkyCoord(ra=tracklet['ra']*u.degree,dec=tracklet['dec']*u.degree, frame="icrs",unit="degree")
    ec_coords=eq_coords.transform_to('geocentrictrueecliptic')  # ecliptic tracklet coordinates
    mean_el = np.mean(ec_coords.lon.value)
    mean_eb = np.mean(ec_coords.lat.value)
# Proper motion in equatorial & ecliptic coords
    fit_mjd = (tracklet['mjd']-mean_mjd)*24                     # tracklet MJDs [hr] (average-subtracted) to fit PM
    fit_ra = (tracklet['ra']-mean_ra)*3600                      # tracklet RAs [''] (average-sub)
    fit_dec = (tracklet['dec']-mean_dec)*3600                   # tracklet Decs [''] (average-sub)
    pmo_ra=cfit(fit_mjd,fit_ra*1000,1,sigma=tracklet['raerr'])
    pmo_dec=cfit(fit_mjd,fit_dec*1000,1,sigma=tracklet['decerr'])
    pm_ra = pmo_ra[0][0]/1000*abs(np.cos(np.deg2rad(mean_dec))) # tracklet ra pm [''/hr?]
    pm_raerr = pmo_ra[1][0]
    pm_dec = pmo_dec[0][0]/1000                                 # tracklet dec pm [''/hr or day?]
    pm_decerr = pmo_dec[1][0]
    fit_el = (np.array(ec_coords.lon.value)-mean_el)*3600       # tracklet ecliptic l's [''] (average-subtracted)
    fit_eb = (np.array(ec_coords.lat.value)-mean_eb)*3600       # tracklet ecliptic b's' [''] (average-sub)
    pmo_el=cfit(fit_mjd,fit_el,1,sigma=None)
    pmo_eb=cfit(fit_mjd,fit_eb,1,sigma=None)
    pm_el = pmo_el[0][0]*abs(np.cos(np.deg2rad(mean_eb)))       # tracklet el pm [''/hr?]
    pm_eb = pmo_eb[0][0]                                        # tracklet eb pm [''/hr or day?]
    redchi = (np.dot(np.transpose(pmo_ra[1]),np.dot(np.linalg.inv(pmo_ra[2]),pmo_ra[1]))) #reduced chi squared
# Magnitudes
    mean_magauto = np.mean(tracklet['mag_auto'])                # mean tracklet magnitude across all bands
    filts=["u","g","r","i","z","Y","VR"]
    filt_array=np.repeat(-99.99,len(filts))                     # for each filter (u g r i z Y VR)
    filt_counter = 0
    for filtr in filts:
        filt_mmts_bool = tracklet['filter']==filtr
        filt_mmts = tracklet[filt_mmts_bool]['mag_auto']
        if len(filt_mmts)>0:
            mean_mag = np.mean(filt_mmts)
            filt_array[filt_counter] = mean_mag
        filt_counter+=1
    return([first_obs_mjd,dt_mjd,tracklet_id_val,mean_ra,mean_dec,mean_el,mean_eb,
            pm_ra,pm_raerr,pm_dec,pm_decerr,pm_el,pm_eb,
            mean_magauto,redchi,nmeas,pix128,pix64,pix32]
           +list(filt_array),pix64,tracklet_id_val,pix32,ring256)


def calc_mpc80(obs,scope_names,scope_codes,cent_letters,cent_numbers):
    '''Calculates necessary info for MPC 80-column format line for observation.
    Arguments:
    ----------
    obs (astropy row)
            a tracklet observation
    scope_names (list)
            names of telescope/observatories used
    scope_codes (list)
            codes for "                          "
    Returns:
    --------
    None; directory "dir" is created if not already there
    '''
#Date of Observation: "YYYY MM DD.ddddd"
    tim=Time(float(obs['mjd']),format="mjd")
    yyyy=str(tim.datetime.year).zfill(4)
    century_marker = cent_letters[cent_numbers==int(yyyy[:2])]
    mm=str(tim.datetime.month).zfill(2)
    dd='{:.5f}'.format(round((tim.datetime.day+(tim.datetime.hour+
        ((tim.datetime.minute+((tim.datetime.second+(tim.datetime.microsecond/1000000))/60))/60))/24),6)).zfill(8)
    obs_date=(century_marker+yyyy+" "+mm+" "+dd)
#Observed RA: "HH MM SS.ddd",
    co=SkyCoord(obs['ra'],obs['dec'],frame="icrs",unit="degree")
    hhr=int(co.ra.hour)
    mmr=int((co.ra.hour-hhr)*60)
    ssr='{:.3f}'.format(round(((co.ra.hour-hhr)*60-mmr)*60,3))
    obs_ra=" ".join([str(hhr).zfill(2),str(mmr).zfill(2),str(ssr).zfill(6)])
#Observed Dec: "sDD MM SS.dd",
    dd=int(co.dec.deg)
    sine=str("+")
    if np.sign(dd)==-1:
        sine=str("-")
    mm=int((co.dec.deg-dd)*60)
    ss='{:.2f}'.format(abs(round((((co.dec.deg-dd)*60-mm)*60),2)))
    obs_dec=''.join([sine,str(abs(dd)).zfill(2)," ",str(abs(mm)).zfill(2)," ",str(ss).zfill(5)])
#Packed Provisional Number:
    #p1="A"
    #for cent in range(0,len(cent_numbers)):
    #   print(int(yyyy[:3]))
    #   if int(yyyy[:3])==cent_numbers[cent]:
    #       p1=cent_letters[cent]
    p_des=str(pix128).zfill(6) #my temporary designation number(the healpix number)
#Observed Magnitude and Band: mm.mmB
    obs_mag=str('{:0.2f}'.format(obs['mag_auto'])).zfill(5)+obs['filter'][0]
#Observatory code:
    instrument=obs['instrument']
    for ins in range(0,len(scope_names)):
            if instrument==scope_names[ins]:
                obs_code=scope_codes[ins]
#Tholen-style uncertainties:
    raerr = str('{:0.3f}'.format(obs['raerr'])).zfill(5)
    decerr = str('{:0.3f}'.format(obs['decerr'])).zfill(5)
    return(obs_date,obs_ra,obs_dec,p_des,obs_mag,obs_code,raerr,decerr)

def make_80col(obs_num,t_num,obs_date,obs_ra,obs_dec,p_des,obs_mag,obs_code,raerr,decerr):
    '''Formats observation info into MPC 80-column line for Find_Orb/d2 use.
    Arguments:
    ----------
    obs_num (int)
            a unique observation number
    t_num (int)
            a unique tracklet number
    the rest is output from calc_mpc80()
    Returns:
    --------
    the observation's MPC 80-column line
    '''
    c01t06=str(obs_num).zfill(2)+"."+str(t_num).zfill(3) #columns "1-5"1-6 observation number, left-padded with 0s
    c07t12="."+str(p_des).zfill(6)   #columns "6-12"7-13 packed provisional number (pix for now)
    c13t14=str(" ")              #columns "13-"14 3X
    c15t32=str(obs_date)+" "     #columns 15-32 Date of observation
    c33t44=str(obs_ra)           #columns 33-44 Observed RA (J2000.0) (deg)
    c45t56=str(obs_dec)          #columns 45-56 Observed Decl. (J2000.0) (deg)
    c57t65=str("         ")      #columns 57-65 9X     Must be blank
    c66t71=str(obs_mag)          #columns 66-71 F5.2,A1   Observed magnitude and band
    c72t77=str("      ")         #columns 72-77 6X      Must be blank
    c78t80=str(obs_code)         #columns 78-80 A3     Observatory code
    c81t86=str(" "+str(raerr))   #columns 81 (blank) and 82-86 = uncertainty in RA (arcsec)
    c87t92=str(" "+str(decerr))  #columns 86 (blank) and 87-92 = uncertainty in Dec (arcsec)
    the_line="".join([c01t06,c07t12,c13t14,c15t32,c33t44,c45t56,c57t65,c66t71,c72t77,c78t80,c81t86,c87t92,str("\n")]) #the 80col line to write!
    return(the_line,c01t06+c07t12)


#---------------------------------------------------------------------------------
# Main Code
#---------------------------------------------------------------------------------
if __name__ == "__main__":

    # --Initiating & Defining--
    subdir = str(sys.argv[1])                         # hgroup subdir#
    drnum = str(sys.argv[2])                          # NSC DR#
    basedir = "/home/x25h971/canfind_dr"+drnum+"/"    # base location of operations
    hpdir = basedir+"hpix/hgroup_"
    concatdir = basedir+"concats/"
    makedir(concatdir)

    # stuff for querying missing mmt info
    qc.set_timeout_request(1800)                      # set datalab query timeout
    mmt_cols = ["measid","objectid","mjd","ra","dec",
                "mag_auto","magerr_auto","filter","raerr","decerr",
                "exposure"] #,"pix128","ring256"]     # previously "ns"
    mmt_dts = ["U","U","float64","float64","float64",
               "float64","float64","U","float64","float64",
               "U"] #,"int","int"]                    # previously "dts"
    filename_exposures = "exp_inst_dr2.fits"          # exposure list for NSC DR2
    file_exposures=Table.read(basedir+"files/"+filename_exposures,format="fits") # exposure table to get instrument column for MPC-80 lines

    # stuff for tracklet info,  tracklet id format will be "subdir.tracklet_#"
    tlet_cols = ["mjd","dmjd","tracklet_id","ra","dec",
                "el","eb","pm_ra","pm_raerr","pm_dec",
                "pm_decerr","pm_el","pm_eb","mag_ave","redchi",
                "nmeas","pix128","pix64","pix32",
                "u","g","r","i","z","Y","VR",]
    tlet_dts = ["float64","float64","U","float64","float64",
                "float64","float64","float64","float64","float64",
                "float64","float64","float64","float64","float64",
                "float64","float64","float64","float64",
                "float64","float64","float64","float64","float64","float64","float64"]

    # stuff for MPC 80-col format info
    scope_codes=["V00","695","W84"]                   # telescope information
    scope_names=["ksb","k4m","c4d"]
    cent_letters=["I","J","K"]                        # century of observation
    cent_numbers=[int(18),int(19),int(20)]

    # set up output files
    cat_mmts_filename = "cf_dr"+drnum+"_hgroup_"+subdir+"_mmts.fits"
    cat_mmts = Table() #names=mmt_cols,dtype=mmt_dts)      # table for tracklet mmts, to be written to output file (1) 
    cat_tracklets_filename = "cf_dr"+drnum+"_hgroup_"+subdir+"_tracklets.fits"
    cat_tracklets = Table(names=tlet_cols,dtype=tlet_dts) # table for tracklet info, to be written to output file (2)
    mpc_filename = "cf_dr"+drnum+"_hgroup_"+subdir+".txt"
    if os.path.exists(concatdir+mpc_filename): os.remove(concatdir+mpc_filename)
    mpc_mmts = open(concatdir+mpc_filename,"a+")          # text file for subdir to store MPC-formatted data, output file (3), previously "file80"
    mpc_mmts.write("123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|12\n")
    mpc_mmts.write("<-ObjDesig->*nnYYYY MM DD.DDDDD HH MM SS.SSSsdd mm ss.ss<blanks >MM.MMBz<ref>COD a.aaa a.aaa\n")

    # --Loop through HP files in subdir--
    for root,dirs,files in os.walk(hpdir+subdir):
        for name in files:
            fname = os.path.join(root,name)
            px = str(name).split(".")[0]
            px = px.split("_")
            pix128 = int(px[1])
            ring256 = -99
            rcount = 0
            if len(px)>2:
                ring256 = px[2] #for dr2
                # check for exsisting tracklets in catalog with same pix128
                rcount = len(cat_tracklets[cat_tracklets['pix128']==pix128])
            if os.stat(os.path.join(root,name)).st_size!=0 and len(fits.open(os.path.join(root,name)))>1:

                # read HP file, get tracklet mmts
                print("reading info from ",fname)
                dat = Table.read(fname)
                hp_mmts = dat[dat['cluster_label']!=-1]       # tracklet mmts for this PIX
                if len(hp_mmts)>0:
                    if ring256==-99: ring256=hp.ang2pix(256,dat[0]['ra'],dat[0]['dec'],lonlat=True)
                    hp_mmts['ring256'] = Column(np.repeat(float(ring256),len(hp_mmts)))
                    hp_cols = np.array(hp_mmts.colnames)

                    # --Query missing columns from nsc.meas--
                    missing = np.invert(np.isin(np.array(mmt_cols),hp_cols)) #boolean array for missing columns
                    missing_cols = np.array(mmt_cols)[missing]
                    missing_dts = np.array(mmt_dts)[missing]
                    if len(missing_cols)>0:
                        missing[0] = True                         # so you always at least have measid
                        qcoltxt = "SELECT "+", ".join(["meas."+(mcol) for mcol in missing_cols])+" FROM nsc_dr"+drnum+".meas as meas JOIN nsc_dr"+drnum+".object as obj on obj.id=meas.objectid WHERE "
                        #print("qcoltxt = "+qcoltxt)
                        if ring256==-99:                      # if this is a full PIX128 file,
                            if int(drnum)==1:                 # if it's NSC DR1, there's no PIX128 column so you gotta get those ring256s!
                                coords128 = hp.pix2ang(128,pix128,lonlat=True)
                                nbs512=hp.get_all_neighbours(512,coords128[0],coords128[1],lonlat=True) # 8 nearest neighbors, nside=512
                                coords512=hp.pix2ang(512,nbs512,lonlat=True)                            # center coordinates for ^^
                                ring256s=np.unique(hp.ang2pix(256,coords512[0],coords512[1],lonlat=True))
                                qtext = qcoltxt+"obj.ring256="+str(ring256s[0])+" or obj.ring256="+str(ring256s[1])+" or obj.ring256="+str(ring256s[2])+" or obj.ring256="+str(ring256s[3])
                            else: qtext = qcoltxt+"obj.PIX="+str(pix128)
                        else: qtext = qcoltxt+"obj.ring256="+str(ring256)                 # if only 1/4 of a full PIX128 file (a ring256)
                        try: dd = qc.query(sql=qtext,fmt="csv")                           # attempt a query!
                        except Exception as e: print("query failed")
                        else: the_query = query_fix(dd,missing_dts)
                        hp_mmts = join(hp_mmts,the_query,keys="measid",join_type="left")  # add queried columns to existing HP tracklet mmt table
                    hp_mmts=join(hp_mmts,file_exposures,keys="exposure",join_type="left") # add exposure info to table
                    cf_output_filename = concatdir+"hgroup128_"+subdir+"/tracklets_"+px[0]+"_"+px[1]+".fits"
                    print("writing ",cf_output_filename)
                    makedir(concatdir+"hgroup128_"+subdir)
                    hp_mmts.write(cf_output_filename,overwrite=True) # write a new fits file with just the tracklet mmts
                    hp_mmts['pix32'] = Column(np.repeat(float(000000),len(hp_mmts)))
                    hp_mmts['pix64'] = Column(np.repeat(float(000000),len(hp_mmts)))
                    hp_mmts['tracklet_id'] = Column(dtype='U100',length=len(hp_mmts))
                    hp_mmts['obs_id'] = Column(dtype='U100',length=len(hp_mmts))
                    if 'pix' in hp_mmts.colnames:
                        hp_mmts['pix128'] = hp_mmts['pix']
                        hp_mmts.remove_column("pix")

                    t_num=0+rcount #counter for tracklet_number (every tracklet in this hpix will have a unique number)
                    hp_mmts.sort("mjd")
                    # --Calculate tracklet info--
                    for trl in np.unique(hp_mmts['cluster_label']):  # for every tracklet,
                        tr=np.where(hp_mmts['cluster_label']==trl)[0]
                        tracklet=hp_mmts[tr]
                        tracklet.sort('mjd')                         # sort in ascending obs_date order
                        #print("cluster ",trl,", len = ",len(tracklet))
                        if len(tracklet)<100 and len(tracklet)>2:
                            print("tracklet = \n",tracklet)
                            tracklet_row,pix64,tid,pix32,r256 = tracklet_calc(tracklet, pix128, ring256,t_num)
                            #print("tracklet number & id = ",t_num,tid)
                            cat_tracklets.add_row(tracklet_row)          # add tracklet info to hgroup table
                            hp_mmts['pix32'][tr] = pix32
                            hp_mmts['pix64'][tr] = pix64
                            hp_mmts['pix128'][tr] = pix128
                            hp_mmts['ring256'][tr] = r256
                            hp_mmts['tracklet_id'][tr] = tid
                            #print("number of measurements = ",len(hp_mmts[tr]))
                            print("tracklet id = ",tid)
                            # --Calculate info for each mmt's MPC 80-col line--
                            obs_num=0 #counter for tracklet observation number
                            for obs in tracklet: #for every observation,
                                obs_80s=calc_mpc80(obs,scope_names,scope_codes,cent_letters,cent_numbers) #obs_80s = [obs_date,obs_ra,obs_dec,obs_mag,obs_code]
                                linne,oid = make_80col(obs_num,t_num,obs_80s[0],obs_80s[1],obs_80s[2],obs_80s[3],obs_80s[4],obs_80s[5],obs_80s[6],obs_80s[7])
                                mpc_mmts.write(linne) # write the line to hgroup MPC80 txt file
                                #print("observation id = ",str(obs_num).zfill(2)+"."+str(t_num).zfill(3)+"."+str(pix128).zfill(6))
                                #print(tr[obs_num]) #,hp_mmts[tr[obs_num]])
                                hp_mmts['obs_id'][tr[obs_num]] = (str(obs_num).zfill(2)+"."+str(t_num).zfill(3)+"."+str(pix128).zfill(6))
                                obs_num+=1
                            t_num+=1
                        #time.sleep(2)

                cat_mmts=vstack([cat_mmts,hp_mmts])                                   # add HP tracklet mmts to hgroup table
    cat_mmts['pix128']=Column(cat_mmts['pix128'],dtype="float64")
    cat_mmts.write(concatdir+cat_mmts_filename,overwrite=True)
    cat_tracklets.write(concatdir+cat_tracklets_filename,overwrite=True)
    mpc_mmts.close()
