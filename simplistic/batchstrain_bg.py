#!/usr/bin/env python

#-----------------------
# Import needed modules
#-----------------------
from __future__ import division
import os
import numpy as np
import pylab as plt
import gwpy
import h5py
import lal
from astropy.coordinates import get_sun
import astropy.time as Time
import astropy
from gwpy.timeseries import TimeSeries
from scipy.signal import butter
from scipy.signal import filtfilt

#------------
# Read macros
#------------
from optparse import OptionParser
import sys
parser = OptionParser()
parser.add_option("--start", dest="starttime", type="int", help="Inputthe GPS ")
parser.add_option("--end", dest="endtime", type="int", help="Input the GPS end ")
(opts, args) = parser.parse_args()
if not opts.starttime:
  print "Error... a '--start' option giving the GPS start time is required"
  sys.exit(1)

starttime = opts.starttime
endtime = opts.endtime
if starttime < 0 or np.isinf(starttime):
  print "Error... the given start time is not valid"
  sys.exit(1)
print starttime
print endtime
#----------------
# Read Timeseries
#----------------
# Read segment instead using gwpy, where starttime and endtime are from the dag file
pathtoinput = "/home/spxha/"
strainH = TimeSeries.read(pathtoinput+'S6framesH1.lcf',channel='H1:LDAS-STRAIN', start=starttime, end=endtime)
strainL = TimeSeries.read(pathtoinput+'S6framesL1.lcf',channel='L1:LDAS-STRAIN', start=starttime, end=endtime)
#-----------
# Xspacing
#-----------
tsH = 2.44140625E-4
tsL = 2.44140625E-4
Xspacing = tsH

#-----------
# Duration
#-----------
gpsStartH = starttime
durationH = endtime - starttime
gpsEndH = endtime

gpsStartL = gpsStartH
durationL = durationH
gpsEndL = gpsEndH

#############################
#---------------------------
# Applying a bandpass filter
#---------------------------
ord = 4
Wn = [100.0/2918.0,300.0/2918.0]
type = 'bandpass'
bL,aL = butter(ord,Wn, btype=type)
bH,aH = butter(ord,Wn, btype=type)
strainL = filtfilt(bL,aL,strainL)
strainH = filtfilt(bH,aH,strainH)

#-----------------------
# Create a time vector
#-----------------------
timeH = np.arange(gpsStartH, gpsEndH, tsH)
timeL = np.arange(gpsStartL, gpsEndL, tsL)
timel=timeL
#-----------------------------------------------
# Create mapping between detector name and index
#-----------------------------------------------
detMap = {'H1': lal.LALDetectorIndexLHODIFF, 'L1':
lal.LALDetectorIndexLLODIFF}

#-------------------------------------
# Get detector structure for H1 and L1
#-------------------------------------
detH1 = lal.CachedDetectors[detMap['H1']]
detL1 = lal.CachedDetectors[detMap['L1']]

#--------------------------------------------------------------
# Set a GPS time (this is just a random value for this example)
#--------------------------------------------------------------
tgps = lal.LIGOTimeGPS(gpsStartH, 0)

#---------------------------------------------------------
# Get right ascension and declination of source in radians
#---------------------------------------------------------
# do this at every 30 seconds
numseg30 = int((endtime-starttime)/30.)
print numseg30
seg30 = gpsStartH + 30*np.linspace(0,numseg30,numseg30)
seg30 = np.append(seg30, endtime)
tdelay = [[0] for _ in range(numseg30)]
for i in range(numseg30):
        if ((timeL[int(i/Xspacing)]>seg30[i])&(timel[int(i/Xspacing)]<seg30[i+1])):
		coordstime=seg30[i]
		coords = get_sun(Time.Time(coordstime,format='gps'))
		tdelay[i] = lal.ArrivalTimeDiff(detH1.location, detL1.location, coords.ra.hour*np.pi/12, coords.dec.hour*np.pi/12, tgps)
	else:
		pass
tdelay = np.repeat(tdelay,int(30/tsH))

# make sure tdelay and timeL are of same length in case integer-ing caused slight inconsistency.
b = np.ones(len(timeL)-len(tdelay))*tdelay[-1]
tdelay = np.append(tdelay,b)

timeL = timeL - tdelay

# at this point timeH and timeL are synchronised
#
#initialise idx, strain30L/H and strainprod
overtsH = int(30/tsH)
# idxH,idxL,strain30L,strain30H = [[[0 for _ in range(overtsH)] for _ in range(numseg30)] for _ in range(4)]
# strainprod = [[0] for _ in range(numseg30)]
# for i in range(numseg30):
# 	idxH = np.where(np.logical_and(seg30[i]<timeH, seg30[i+1]>timeH))
# 	idxL = np.where(np.logical_and(seg30[i]<timeL, seg30[i+1]>timeL))
# 	strain30L = strainL[idxL]
# 	strain30H = strainH[idxH]
# 	strainprod[i] = np.dot(strain30H,strain30L)

#---------------
# Get Background
#---------------
background_intervals = np.linspace(1,50,50)*600
ra_background, background_tdelay, dec_background, coords_background = [[[0] for _ in range(50)] for _ in range(4)]
background_intervals = background_intervals.astype(int)
for j in range(len(background_intervals)):
	coords_background[j] = get_sun(Time.Time(gpsStartH-background_intervals[j],format='gps'))
	ra_background[j]     = coords_background[j].ra.hour  * np.pi/12
	dec_background[j]    = coords_background[j].dec.hour * np.pi/12
	background_tdelay[j] = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra_background[j], dec_background[j], tgps)
# initialise backtimeH and backtimeL
backtimeH, backtimeL = [[[0] for _ in range(50)] for _ in range(2)]
for i in range(50):
        backtimeH[i] = timeH
        backtimeL[i] = timeL - background_tdelay[i]

bg_idxH,bg_idxL,bg_strain30L,bg_strain30H = [[[0 for _ in range(overtsH)] for _ in range(numseg30)] for _ in range(4)]
bg_strainprod = [[[0] for _ in range(numseg30)] for _ in range(len(background_intervals))]
for j in range(len(background_intervals)):
	for i in range(numseg30):
		bg_idxH = np.where(np.logical_and(seg30[i]<backtimeH[j], seg30[i+1]>backtimeH[j]))
		bg_idxL = np.where(np.logical_and(seg30[i]<backtimeL[j], seg30[i+1]>backtimeL[j]))
		bg_strain30L = strainL[bg_idxL]
		bg_strain30H = strainH[bg_idxH]
		while len(bg_strain30H)>len(bg_strain30L):
                  bg_strain30H = np.delete(bg_strain30H,-1)
                while len(bg_strain30H)<len(bg_strain30L):
                  bg_strain30L = np.delete(bg_strain30L,-1)
                bg_strainprod[j][i] = np.dot(bg_strain30H,bg_strain30L)


#---------------------------
# Write output to hdf5 file
#---------------------------
# outputf = h5py.File(str(starttime)+'.hdf5','w')
# siggroup = outputf.create_group('Signal')
# timesetH = siggroup.create_dataset('TimeSeries',data=timeH)
# for i in range(len(strainprod)):
# 	siggroup.create_dataset('Strain'+str(i),data=strainprod[i])

outdir='/scratch/spxha/'+str(starttime)+'_'+str(endtime)
bg_group = [[0] for _ in range(len(background_intervals))]
if os.path.exists(outdir)==True:
        print "True"
	outputf = h5py.File(outdir+'/'+str(starttime)+'.hdf5','r+')
	for j in range(len(background_intervals)):
		bg_group[j] = outputf.create_group('Background'+str(j))
		for i in range(numseg30):
			print i
			bg_group[j].create_dataset('Strain'+str(i),data=bg_strainprod[j][i])
else:
	pass

