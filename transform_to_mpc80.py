
#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# transform_to_mpc80.py is a script that reads fits files in a subdir (hgroup_<subdir#>) with CANFind output,
# grabs the tracklets, and records them in the Minor Planet Center's 80-column format,
# explained in no detail below (find more at https://www.minorplanetcenter.net/iau/info/OpticalObs.html),
# then writes the new tracklet table to a text file.  This is the file type that digest2 and MPC require.
# It also queries any necessary information? We'll see about that. The tracklet measurements are also
# written to a FITS file, to later be concatenated into the tracklet catalog by cf_tracklet_cat_creator.py

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
#    14            A1     Note 1
#    15            A1     Note 2
#    16 - 32              Date of observation
#    33 - 44              Observed RA (J2000.0)
#    45 - 56              Observed Dec (J2000.0)
#    57 - 65       9X     Must be blank
#    66 - 71    F5.2,A1   Observed magnitude and band
#    72 - 77       X      Must be blank
#    78 - 80       A3     Observatory code

# Example Row:
#
#
#

#---------------------------------------------------------------------------------
# Imports
#---------------------------------------------------------------------------------
import matplotlib
import sys
import os
import numpy as np
from astropy.table import Table,Column,join,vstack
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.time import Time
from dl import queryClient as qc
import healpy as hp
#---------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------
def query_fix(quer,dtps):
	bc=[]
	bdd=quer.splitlines()
	for ag in bdd:
		bc.append(ag.split(","))
	return(Table(np.array(bc[1:]),names=np.array(bc[0]),dtype=dtps))

def makedir(dir):
	'''makes a directory with name "dir" if it does not exist'''
	if not os.path.exists(dir):
		os.mkdir(dir)
#---------------------------------------------------------------------------------
# Main Code
#---------------------------------------------------------------------------------
if __name__ == "__main__":

	scope_codes=["V00","695","W84"]
	scope_names=["ksb","k4m","c4d"]
	cent_letters=["I","J","K"]
	cent_numbers=[int(18),int(19),int(20)]

	ns=["mjd","ra","dec","measid","objectid","mag_auto","filter","raerr","decerr","exposure"]
	dts=["float64","float64","float64","U","U","float64","U","float64","float64","U"]
	cat_mmts=Table(names=ns,dtype=dts)
	qc.set_timeout_request(1800)

	subdir=sys.argv[1] #the subdirectory number containing the files to combine (hgroup_subdir)
	makedir("concat_files/hgroup_%s"%subdir)
	t_num=0 #a counter for tracklet_number (every tracklet in this subdir hgroup will have a unique number)
	file80=open("canfind_tracklet_obs_%s.txt" % subdir,"a+") #write the text file with filenumber to store data
	file_exposures=Table.read("exp_inst.fits",format="fits") #NSC DR1 exposure table 
	for root,dirs,files in os.walk("hgroup_%s"%str(subdir)):#tracklet mmts from files in subdir will be witten to text file. 
		for name in files:
			if os.stat(os.path.join(root,name)).st_size!=0 and len(fits.open(os.path.join(root,name)))>1:
				hdul=fits.open(os.path.join(root,name))
				#print("name=",os.path.join(root,name))
				dat=Table(hdul[1].data)
				#dat.sort("mjd")
				hdul.close()
				#dat=Table.read(os.path.join(root,name),format="fits") #read a hpix FITS file
				bool_out=dat['cluster_label']!=-1
				t_out=dat[bool_out]
				#if there are actually tracklets in this healpix,
				#print(len(t_out))
				pix=str(name)[8:-5]
				if len(t_out)>0 and (not (str(pix).zfill(6)) in file80.read()):
				#if len(t_out)>0 and not (str(pix).zfill(6)) in file80.read() and not os.path.exists("concat_files/cf_tracklets_hp_%s.fits" % str(pix)):  
					#Do a query to get the errors and instrument information
					RA=hp.pix2ang(128,int(pix),lonlat=True)[0]
					DEC=hp.pix2ang(128,int(pix),lonlat=True)[1]
					nbs=hp.get_all_neighbours(512,RA,DEC,lonlat=True) #get the 8 nearest neighbors to the cooordinates for nside=512
					coords=hp.pix2ang(512,nbs,lonlat=True) #get the center coordinates for the8 nside=512  neighbors
					fpix=np.unique(hp.ang2pix(256,coords[0],coords[1],lonlat=True))
					qnames=["measid","objectid","raerr","decerr","exposure","filter","mag_auto"]
					qdtypes=['U','U','f8','f8','U','U','f8']
					qtext="".join(["SELECT meas.measid,meas.objectid,meas.raerr,meas.decerr,meas.exposure,meas.filter,meas.mag_auto FROM nsc_dr1.meas as meas JOIN nsc_dr1.object as obj on obj.id=meas.objectid WHERE obj.ring256=",str(fpix[0])," or obj.ring256=",str(fpix[1])," or obj.ring256=",str(fpix[2])," or obj.ring256=",str(fpix[3])]) 
					try:
						dd=qc.query(sql=qtext,fmt="csv")
					except Exception as e:
						if "Query timeout" in e.message: #if the query timed out,
							the_query=Table(names=qnames,dtype=qdtypes)
							for iq in [0,1,2,3]:
								qqtext="".join(["SELECT meas.measid,meas.objectid,meas.raerr,meas.decerr,meas.exposure,meas.filter,meas.mag_autoFROM nsc_dr1.meas as meas JOIN nsc_dr1.object as obj on obj.id=meas.objectid WHERE obj.ring256=",str(fpix[iq])]) #split query into 4, try again
								ddq=qc.query(sql=qtext,fmt="csv")
								the_query_temp=query_fix(ddq,qdtypes)
								the_query=vstack([the_query,the_query_temp])
						else:
							print("query failed")
					else:
						the_query=query_fix(dd,qdtypes)
					#add the queried information to the tracklet table! we'll write it to a fits file.
					if ("objectid" in t_out.colnames):
						t_out.remove_column("objectid")
					if ("mag_auto" in t_out.colnames):
						t_out.remove_column("mag_auto")
					if ("filter" in t_out.colnames):
						t_out.remove_column("filter")
					t_out=join(t_out,the_query,keys="measid",join_type="left")
					t_out=join(t_out,file_exposures,keys="exposure",join_type="left")
					#write a new fits file with just the tracklet mmts
					t_out.write("concat_files/hgroup_%s/cf_tracklets_hp_%s.fits" % (subdir,pix),format="fits",overwrite=True)
					cat_mmts=vstack([cat_mmts,t_out])
				#information with added columns for the healpix
				#for every tracklet:
					for trl in np.unique(t_out['cluster_label']):
						tr=t_out['cluster_label']==trl
						tracklet=t_out[tr]
						obs_num=0 #a counter for observation number within the traklet
					#for every observation in that tracklet:
						for obs in tracklet:
						#Necessary Calculations
						#Date of Observation: "YYYY MM DD.dddddd"
							tim=Time(float(obs['mjd']),format="mjd")
							yyyy=str(tim.datetime.year).zfill(4)
							mm=str(tim.datetime.month).zfill(2)
							dd='{:.6f}'.format(round((tim.datetime.day+(tim.datetime.hour+
								((tim.datetime.minute+((tim.datetime.second+
									(tim.datetime.microsecond/1000000)
									)/60))/60)
								)/24),6))
							obs_date=" ".join([yyyy,mm,dd.zfill(9)])
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
							obs_dec="".join([sine,str(abs(dd))," ",str(abs(mm)).zfill(2)," ",str(ss).zfill(5)])
						#Packed Provisional Number:
							#p1="A"
							#for cent in range(0,len(cent_numbers)):
							#	print(int(yyyy[:3]))
							#	if int(yyyy[:3])==cent_numbers[cent]:
							#		p1=cent_letters[cent]
							#p_des=str(pix).zfill(6) #my temporary designation number(the healpix number)
						#Observed Magnitude and Band: mm.dd,I
							obs_mag="".join([str(round(obs['mag_auto'],2)).zfill(5),obs['filter'][0]])
						#Observatory code:
							instrument=obs['instrument']
							for ins in range(0,len(scope_names)):
	    							if instrument==scope_names[ins]: 
	        							obs_code=scope_codes[ins]
						#Create the line in 80-column format to write to the txt file
							c01t05=str(obs_num).zfill(6) #columns 1-5 observation number, left-padded with 0s
							c06t12=str(t_num).zfill(7)   #columns 6-12 packed provisional number (pix for now)
							c13t15=str("   ")            #columns 13-15 3X
							c16t32=str(obs_date)         #columns 16-32 Date of observation
							c33t44=str(obs_ra)           #columns 33-44 Observed RA (J2000.0)
							c45t56=str(obs_dec)          #columns 45-56 Observed Decl. (J2000.0)
							c57t65=str("         ")      #columns 57-65 9X     Must be blank
							c66t71=str(obs_mag)          #columns 66-71 F5.2,A1   Observed magnitude and band
							c72t77=str("      ")         #columns 72-77 6X      Must be blank
							c78t80=str(obs_code)         #columns 78-80 A3     Observatory code
							linne="".join([c01t05,c06t12,c13t15,c16t32,c33t44,c45t56,c57t65,c66t71,c72t77,c78t80,str("\n")]) #the 80col line to write!
							file80.write(linne)          #write the line
							obs_num+=1
						t_num+=1
	file80.close()
	cat_mmts.write("concat_files/cat_files/canfind_mmts_hgroup_%s.fits"%subdir,overwrite=True)
	os.replace("canfind_tracklet_obs_%s.txt" % subdir, "concat_files/txt_files/canfind_tracklet_obs_%s.txt" % subdir)



