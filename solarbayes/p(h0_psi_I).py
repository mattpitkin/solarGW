#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy, gwpy, h5py, lal
from astropy.coordinates import get_sun
import astropy.time as Time
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res


#-------- First importing some stuff from earlier work --------#
starttime=931219808
endtime = 971614889

gpsStartH = starttime
durationH = endtime - starttime
gpsEndH = endtime

gpsStartL = gpsStartH
durationL = durationH
gpsEndL = gpsEndH

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
# do this every 30 seconds
numseg30 = int((endtime-starttime)/30.)
seg30 = gpsStartH + 30*np.linspace(1,numseg30,numseg30)
tdelay = [[0] for _ in range(numseg30)]
for i in range(numseg30-1):
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

# H1 and L1 are now in sync.

#------- Preparing for fc and fp ------#
numseg = np.int((dtime)/600)
segs = np.linspace(1,numseg,numseg)*600
segs = segs + starttime - 600
ra,dec,fp,fc = [[0 for _ in range(numseg)] for _ in range(4)]
for i in range(numseg):
	coordstime = segs[i]
	coords = get_sun(Time.Time(coordstime,format='gps'))
	ra[i] = coords.ra.hour*np.pi/12
	dec[i] = coords.dec.hour*np.pi/12
psi_array = np.linspace(0,np.pi,10)

# X is H1 and Y is L1

sigmaX = np.std(StrainH1) #just std of strainH1
sigmaY = np.std(StrainL1) # std of strainL!
sigmaA = 10.0
h0 = 3*sigmaA

Xspacing = 2.44140625E-4

###########
# Next I want to get

# for i in range(time/Xspacing):
# use function directly for fc and fp

# for psi in psi_array = np.linspace(0,np.pi,10)

# Finding probability distribution
d = np.array([dX, dY])
M = h0*np.array([[FpX, FpY], [FcX, FcY]])
C = np.array([[sigmaX**2, 0.], [0., sigmaY**2]])
invC = np.array([[(1./sigmaX**2), 0.], [0., (1/sigmaY**2)]])
detC = sigmaX**2 * sigmaY**2
invSigma0 = np.array([[(1./sigmaA**2) 0.], [0., (1./sigmaA**2)]])
detSigma0 = sigmaA**4
invSigma = np.dot(M.T, np.dot(invC, M)) + invSigma0
Sigma = np.linalg.inv(invSigma)
detSigma = np.linalg.det(Sigma)
chi = np.dot(Sigma, np.dot(M.T, np.dot(invC, d)))
p = 0.5*np.log(detSigma) - 0.5*log(16.*np.pi**4*detSigma0*detC) -
    0.5*(np.vdot(d.T, np.dot(invC, d)) + np.vdot(chi.T, np.dot(invSigma, chi)))

