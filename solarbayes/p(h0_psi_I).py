#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy, gwpy, h5py, lal
from astropy.coordinates import get_sun
import astropy.time as Time
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res

#-------- Filtering and timeshifting detectors --------#
starttime = 931219808
endtime   = 971614889
gpsStartH = starttime
durationH = endtime - starttime
gpsEndH   = endtime
gpsStartL = gpsStartH
durationL = durationH
gpsEndL   = gpsEndH
#---------------------------
# Applying a bandpass filter
#---------------------------
ord = 4
Wn = [100.0/2918.0,150.0/2918.0]
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
#---------------
# Set a GPS time
#---------------
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

# H1 and L1 are now in sync and filtered between 100 and 150 Hz.

#------- Defining some stuff for p ------#
numseg = int((durationH)/600)
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

sigmaX = np.std(strainH) # just std of strainH1
sigmaY = np.std(strainL) # std of strainL!
sigmaA = 10.0
h0 = np.lin(0,3*sigmaA,10)
invSigma0 = np.array([[(1./sigmaA**2) 0.], [0., (1./sigmaA**2)]])
detSigma0 = sigmaA**4
dX = strainH
dY = strainY

###########
# Finding probability distribution
# for psi in psi_array
d, invC, detC, invSigma, detSigma, chi = [[0 for _ in range(int(durationH/Xspacing))] for _ in range(6)]
for i in range(int(durationH/Xspacing)):
	psi = np.pi/2.0
	FpX, FcX = ant_res(gpsTime[i], ra[i], dec[i], psi, 'H1')
	FpY, FcY = ant_res(gpsTime[i], ra[i], dec[i], psi, 'L1')
	d[i] = np.array([dX[i], dY][i])
	M = h0*np.array([[FpX, FpY], [FcX, FcY]])
	C = np.array([[sigmaX**2, 0.], [0., sigmaY**2]])
	invC[i] = np.array([[(1./sigmaX**2), 0.], [0., (1/sigmaY**2)]])
	detC[i] = sigmaX**2 * sigmaY**2
	invSigma[i] = np.dot(M.T, np.dot(invC[i], M)) + invSigma0
	Sigma = np.linalg.inv(invSigma[i])
	detSigma[i] = np.linalg.det(Sigma)
	chi[i] = np.dot(Sigma, np.dot(M.T, np.dot(invC, d)))
p = 0.5*np.log(detSigma) - 0.5*log(16.*np.pi**4*detSigma0*detC) -
    0.5*(np.vdot(d.T, np.dot(invC, d)) + np.vdot(chi.T, np.dot(invSigma, chi)))

#------ plot the probability distribution
fname = 'probdist.pdf'
with PdfPages(fname) as pdf:
	fig1 = plt.figure()
	plt.plot(h0,p)
	plt.title('Probability Distribution')
	pdf.savefig(fig1)
	plt.close()
print 'Plot saved as', fname, 'Did we find GWs?!'
