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

####################################
# Finding probability distribution #
#------- Defining some stuff for p ------#
numseg = int((durationH)/600)
print numseg
segs = np.linspace(1,numseg,numseg)*600
segs = segs + starttime - 600
ra,dec,fp,fc = [[0 for _ in range(numseg)] for _ in range(4)]
for i in range(numseg+1):
	print i
	coordstime = segs[i]
	coords = get_sun(Time.Time(coordstime,format='gps'))
	ra[i] = coords.ra.hour*np.pi/12
	dec[i] = coords.dec.hour*np.pi/12
psi_array = np.linspace(0,np.pi,10)

# X is H1 and Y is L1

sigmaA = 10.0
h0_array = np.linspace(0,3*sigmaX,25)
invSigma0 = np.array([[(1./sigmaA**2), 0.], [0., (1./sigmaA**2)]])
detSigma0 = sigmaA**4
dX = strainH
dY = strainL
FcX, FpX, FcY, FpY = [[0 for _ in range(int(duration/Xspacing))] for _ in range(4)]
for i in range(int(duration/Xspacing)):
	FpX0[i], FcX0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'H1')
	FpY0[i], FcY0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'L1')


p = [[[0 for _ in range(int(durationH/Xspacing))] for _ in range(len(psi_array))] for _ in range(len(h0_array))]
for j in range(len(h0_array)):
	for k in range(len(psi_array)):
		for i in range(int(durationH/Xspacing)):
			print i
			cos2pi = np.cos(2*psi_array[i])
			sin2pi = np.sin(2*psi_array[i])
			FpX = FpX0[i]*cos2pi + FcX0[i]*sin2pi
			FcX = FcX0[i]*cos2pi - FpX0[i]*sin2pi
			FpY = FpY0[i]*cos2pi + FcY0[i]*sin2pi
			FcY = FcY0[i]*cos2pi - FpY0[i]*sin2pi
			int0 = i - int(30/Xspacing)
			int1 = i + int(30/Xspacing)
			sigmaX = np.std(strainH[int0:int1])
			sigmaY = np.std(strainL[int0:int1])
			d = np.array([dX, dY])
			d.shape = (2,1)
			M = h0_array[j]*np.array([[FpX, FpY], [FcX, FcY]])
			C = np.array([[sigmaX**2, 0.], [0., sigmaY**2]])
			invC = np.array([[(1./sigmaX**2), 0.], [0., (1/sigmaY**2)]])
			detC = sigmaX**2 * sigmaY**2
			invSigma = np.dot(M.T, np.dot(invC, M)) + invSigma0
			Sigma = np.linalg.inv(invSigma)
			detSigma = np.linalg.det(Sigma)
			chi = np.dot(Sigma, np.dot(M.T, np.dot(invC, d)))
			p[k][j][i] = 0.5*np.log(detSigma) - 0.5*log(16.*np.pi**4*detSigma0*detC) -  0.5*(np.vdot(d.T, np.dot(invC, d)) + np.vdot(chi.T, np.dot(invSigma, chi)))

#------ plot the probability distribution
fname = 'probdist.pdf'
with PdfPages(fname) as pdf:
	fig1 = plt.figure()
	plt.plot(h0,p)
	plt.title('Probability Distribution')
	pdf.savefig(fig1)
	plt.close()
print 'Plot saved as', fname, 'Did we find GWs?!'
