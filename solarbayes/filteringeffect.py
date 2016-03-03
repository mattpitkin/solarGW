# Making a simple plot to compare simple filtering to notch-filtering
print 'Importing way too many packages. I really do not need most of them but I am too lazy to check which'
import numpy as np
import astropy, gwpy, h5py, lal
from astropy.coordinates import get_sun
import astropy.time as Time
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res
from scipy.misc import logsumexp
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer, AdaptiveETA, AbsoluteETA, AdaptiveTransferSpeed
from notchfilt import get_filter_coefs, filter_data, iir_bandstops

print 'Reading data'
starttime = 969062862.0
endtime = 969063629.0
gpsStartH = starttime
durationH = endtime - starttime
gpsEndH   = endtime
gpsStartL = gpsStartH
durationL = durationH
gpsEndL   = gpsEndH
Xspacing = 2.44140625E-4
gpsTime = np.linspace(starttime,endtime,int(durationH/Xspacing)*4)
pathtoinput = "/home/spxha/"
strainH = TimeSeries.read(pathtoinput+'S6framesH1.lcf',channel='H1:LDAS-STRAIN', start=starttime, end=endtime)
strainL = TimeSeries.read(pathtoinput+'S6framesL1.lcf',channel='L1:LDAS-STRAIN', start=starttime, end=endtime)
print len(strainH)
num_points = int(durationH/Xspacing)
timeH = np.arange(starttime, endtime, Xspacing)
timeL = np.arange(starttime, endtime, Xspacing)
print len(timeH)
print 'Filtering data'
ord = 4
Wn = [100.0/2918.0,150.0/2918.0]
type = 'bandpass'
bL,aL = butter(ord,Wn, btype=type)
bH,aH = butter(ord,Wn, btype=type)
strainL1 = filtfilt(bL,aL,strainL)
strainH1 = filtfilt(bH,aH,strainH)
print 'Standard deviations before adding the notch-filters', np.std(strainL1), np.std(strainH1)


coefsL = get_filter_coefs('L1')
coefsH = get_filter_coefs('H1')
strainL2 = filter_data(strainL,coefsL)
strainH2 = filter_data(strainH,coefsH)
print 'Standard deviations after adding the notch-filters', np.std(strainL2), np.std(strainH2)
print len(timeL),len(strainL),len(gpsTime)
print 'Plotting'
plt.figure(1)
plt.plot(gpsTime,strainH1,'-')
plt.title('Bandpass filtered data of H1 detector for '+str(durationH/60)+' minutes')
plt.xlabel('TimeSeries')
plt.ylabel('Strain')
plt.savefig('1.png')
plt.show(1)

plt.figure(2)
plt.plot(gpsTime,strainH2,'-')
plt.xlabel('TimeSeries')
plt.ylabel('Strain')
plt.title('Bandpass + notch filtered data of H1 detector for '+str(durationH/60)+' minutes')
plt.savefig('2.png')
plt.show(2)

plt.figure(3)
plt.plot(gpsTime,strainL1,'-')
plt.title('Bandpass filtered data of L1 detector for '+str(durationH/60)+' minutes')
plt.xlabel('TimeSeries')
plt.ylabel('Strain')
plt.savefig('3.png')
plt.show(3)

plt.figure(4)
plt.plot(gpsTime,strainL2,'-')
plt.xlabel('TimeSeries')
plt.ylabel('Strain')
plt.title('Bandpass + notch filtered data of L1 detector for '+str(durationH/60)+' minutes')
plt.savefig('1.png')
plt.show(4)
