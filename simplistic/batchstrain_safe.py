#-----------------------
# Import needed modules
#-----------------------
from __future__ import division
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
starttime=931219808
endtime=starttime+70
#---------------
# Read Timeseries
#---------------
# Read segment instead using gwpy, where starttime and endtime are from the dag file
pathtoinput = "/home/spxha/“
strainH = TimeSeries.read(pathtoinput+'S6framesH1.lcf',channel='H1:LDAS-STRAIN', start=starttime, end=endtime)
strainL = TimeSeries.read(pathtoinput+'S6framesL1.lcf',channel='L1:LDAS-STRAIN', start=starttime, end=endtime)
#---------------------
# Xspacing
#---------------------
tsH = 2.44140625E-4
tsL = 2.44140625E-4
Xspacing = tsH

#-----------------------
# Duration
#-----------------------
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

#---------------------------
# Create a time vector
#---------------------------
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
seg30 = gpsStartH + 30*np.linspace(1,numseg30,numseg30)
tdelay = [[0] for _ in range(numseg30)]
for i in range(numseg30):
        if ((timel[int(i/Xspacing)]>seg30[i])&(timel[int(i/Xspacing)]<seg30[i+1])):
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
#initialise idx, strain30L/H and strainprod
overtsH = int(30/tsH)
idxH,idxL,strain30L,strain30H = [[[0 for _ in range(overtsH)] for _ in range(numseg30)] for _ in range(4)]
strainprod = [[0] for _ in range(numseg30)]]
for i in range(numseg30):
	idxH[i] = np.where(np.logical_and(seg30[i]<timeH, seg30[i+1]>timeH))
	print idxH
	idxL[i] = np.where(np.logical_and(seg30[i]<timeL, seg30[i+1]>timeL))
	strain30L[i] = strainL[idxL[i]]
	strain30H[i] = strainH[idxH[i]]
	print strain30L
	strainprod[i] = np.dot(strain30H[i],strain30L[i])
#---------------
# Get Background
#---------------
background_intervals = np.linspace(1,120,120)*600
ra_background, background_tdelay, dec_background, coords_background = [[[0] for _ in range(120)] for _ in range(4)]
background_intervals = background_intervals.astype(int)
for j in range(len(background_intervals)):
	coords_background[j]=get_sun(Time.Time(gpsStartH-background_intervals[j],format='gps'))
	ra_background[j]  = coords_background[j].ra.hour  * np.pi/12
	dec_background[j] = coords_background[j].dec.hour * np.pi/12
	background_tdelay[j] = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra_background[j], dec_background[j], tgps)

#-------------------------------------------------------
# Write output to hdf5 file
#-------------------------------------------------------
outputf = h5py.File(str(starttime)+'.hdf5','w')
siggroup = outputf.create_group('Signal')
backgroundgroup = outputf.create_group('Background')
timesetH = siggroup.create_dataset('TimeSeries',data=timeH)
# for i in range(len(strainprod)):
siggroup.create_dataset('Strain',data=strainprod)
backgroundset = backgroundgroup.create_dataset('tdelay',data=background_tdelay)
