#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy, gwpy, h5py, lal
from astropy.coordinates import get_sun
import astropy.time as Time
from scipy.signal import butter
from scipy.signal import filtfilt
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res

################################################################
#-------- Importing, filtering and timeshifting data ----------#
starttime = 969062862
endtime   = 969063629
gpsStartH = starttime
durationH = endtime - starttime
gpsEndH   = endtime
gpsStartL = gpsStartH
durationL = durationH
gpsEndL   = gpsEndH
Xspacing = 2.44140625E-4
gpsTime = np.linspace(starttime,endtime,int(1/Xspacing))
pathtoinput = "/home/spxha/"
strainH = TimeSeries.read(pathtoinput+'S6framesH1.lcf',channel='H1:LDAS-STRAIN', start=starttime, end=endtime)
strainL = TimeSeries.read(pathtoinput+'S6framesL1.lcf',channel='L1:LDAS-STRAIN', start=starttime, end=endtime)

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

timeH = np.arange(gpsStartH, gpsEndH, Xspacing)
timeL = np.arange(gpsStartL, gpsEndL, Xspacing)
timel=timeL
detMap = {'H1': lal.LALDetectorIndexLHODIFF, 'L1':
lal.LALDetectorIndexLLODIFF}
detH1 = lal.CachedDetectors[detMap['H1']]
detL1 = lal.CachedDetectors[detMap['L1']]
tgps = lal.LIGOTimeGPS(gpsStartH, 0)

#---------- Get right ascension and declination of source in radians ----------#
numseg30 = int((endtime-starttime)/30.)
seg30 = gpsStartH + 30*np.linspace(1,numseg30,numseg30) # 30 second update rate
tdelay = [[0] for _ in range(numseg30)]
for i in range(numseg30-1):
	print i
        if ((timel[int(i/Xspacing)]>seg30[i])&(timel[int(i/Xspacing)]<seg30[i+1])):
		coordstime=seg30[i]
		coords = get_sun(Time.Time(coordstime,format='gps'))
		tdelay[i] = lal.ArrivalTimeDiff(detH1.location, detL1.location, coords.ra.hour*np.pi/12, coords.dec.hour*np.pi/12, tgps)
	else:
		pass
tdelay = np.repeat(tdelay,int(30/Xspacing))

# make sure tdelay and timeL are of same length in case integer-ing caused slight inconsistency.

b = np.ones(len(timeL)-len(tdelay))*tdelay[-1]
tdelay = np.append(tdelay,b)

timeL = timeL - tdelay

# H1 and L1 are now in sync and filtered between 100 and 150 Hz.

#------- Defining some stuff for p ------#
numseg = int((durationH)/600)
print numseg
segs = np.linspace(1,numseg,numseg)*600
segs = segs + starttime - 600
ra,dec,fp,fc = [[0 for _ in range(numseg)] for _ in range(4)]
for i in range(numseg):
	print i
	coordstime = segs[i]
	coords = get_sun(Time.Time(coordstime,format='gps'))
	ra[i] = coords.ra.hour*np.pi/12
	dec[i] = coords.dec.hour*np.pi/12
psi_array = np.linspace(0,np.pi,10)

# X is H1 and Y is L1

sigmaX = np.std(strainH)
sigmaY = np.std(strainL)
sigmaA = 10.0
# h0 = np.linspace(0,3*sigmaA,10)
h0 = 3*sigmaA
invSigma0 = np.array([[(1./sigmaA**2), 0.], [0., (1./sigmaA**2)]])
detSigma0 = sigmaA**4
dX = strainH
dY = strainL

###########
# Finding probability distribution
# for psi in psi_array
detC, detSigma, chi = [[0 for _ in range(int(durationH/Xspacing))] for _ in range(3)]
d = [[0 for _ in range(2)] for _ in range(int(durationH/Xspacing))]
invC, invSigma = [[[[0 for _ in range(2)] for _ in range(2)] for _ in range(int(durationH/Xspacing))] for _ in range(2)]
psi = np.pi/2.0
for i in range(int(durationH/Xspacing)):
	print i
	FpX, FcX = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], psi, 'H1')
	FpY, FcY = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], psi, 'L1')
	FpX = FpX[0]
	FcX = FcX[0]
	FpY = FpY[0]
	FcY = FcY[0]
	print 'FpX is ...',FpX
	d[i] = np.array([dX[i], dY[i]])
	d[i].shape = (2,1)
	M = h0*np.array([[FpX, FpY], [FcX, FcY]])
	C = np.array([[sigmaX**2, 0.], [0., sigmaY**2]])
	invC[i] = np.array([[(1./sigmaX**2), 0.], [0., (1/sigmaY**2)]])
	detC[i] = sigmaX**2 * sigmaY**2
	invSigma[i] = np.dot(M.T, np.dot(invC[i], M)) + invSigma0
	print 'invSigma[i] is ...',invSigma[i]
	Sigma = np.linalg.inv(invSigma[i])
	detSigma[i] = np.linalg.det(Sigma)
	chi[i] = np.dot(Sigma, np.dot(M.T, np.dot(invC[i], d[i])))
p = 0.5*np.log(detSigma) - 0.5*log(16.*np.pi**4*detSigma0*detC) -  0.5*(np.vdot(d.T, np.dot(invC, d)) + np.vdot(chi.T, np.dot(invSigma, chi)))

#------ plot the probability distribution
fname = 'probdist.pdf'
with PdfPages(fname) as pdf:
	fig1 = plt.figure()
	plt.plot(h0,p)
	plt.title('Probability Distribution')
	pdf.savefig(fig1)
	plt.close()
print 'Plot saved as', fname, 'Did we find GWs?!'
