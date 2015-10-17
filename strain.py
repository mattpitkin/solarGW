#-----------------------
# Import needed modules
#-----------------------
import numpy as np
import pylab as plt
import gwpy
import h5py
import lal
from astropy.coordinates import get_sun
import astropy.time as Time
import astropy
from gwpy.timeseries import TimeSeries

#---------------
# Open the File
#---------------
fileNameH = 'H-H1_LOSC_4_V1-967966720-4096.hdf5'
dataFileH = h5py.File(fileNameH, 'r')

#-------------------
# Explore the file
#-------------------
for key in dataFileH.keys():
    print key

#---------------------
# Read in strain data
#---------------------
strainH = dataFileH['strain']['Strain'].value
tsH = dataFileH['strain']['Strain'].attrs['Xspacing']

#-----------------------
# Read in some meta data
#-----------------------
metaKeysH = dataFileH['meta'].keys()
metaH = dataFileH['meta']
print("\n\n")
for key in metaKeysH:
    print(key, metaH[key].value)

gpsStartH = metaH['GPSstart'].value
durationH = metaH['Duration'].value
gpsEndH   = gpsStartH + durationH
detectorH = metaH['Detector']
dataFileH.close()


#--------------------------#
# Redo this for L detector #
#--------------------------#
fileNameL = 'L-L1_LOSC_4_V1-967966720-4096.hdf5'
dataFileL = h5py.File(fileNameL, 'r')
strainL = dataFileL['strain']['Strain'].value
tsL = dataFileL['strain']['Strain'].attrs['Xspacing']
metaKeysL = dataFileL['meta'].keys()
metaL = dataFileL['meta']
gpsStartL = metaL['GPSstart'].value
durationL = metaL['Duration'].value
gpsEndL   = gpsStartL + durationL
detectorL = metaL['Detector']
dataFileL.close()

#############################
#---------------------------
# Applying a highpass filter
#---------------------------
strainL = strainL.highpass(100)
strainH = strainH.highpass(100)

#---------------------------
# Applying a lowpass filter
#---------------------------
strainL = strainL.lowpass(300)
strainH = strainH.lowpass(300)


#---------------------------
# Create a time vector
#---------------------------
timeH = np.arange(gpsStartH, gpsEndH, tsH)
timeL = np.arange(gpsStartL, gpsEndL, tsL)

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
coords = get_sun(Time(967966720.0,format='gps'))
ra  = coords.ra.hour  * np.pi/12
dec = coords.dec.hour * np.pi/12

#---------------
# Get time delay
#---------------
tdelay = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra, dec, tgps)

#----------------------
# Plot the time series
#----------------------
numSamples = 10000
plt.plot(timeL[0:numSamples]-tdelay, strainL[0:numSamples],label='L1 data')
plt.plot(timeH[0:numSamples], strainH[0:numSamples], label='H1 data')
plt.xlabel('GPS Time (s)')
plt.ylabel('Strain')
plt.title('Strain data on 2010-09-08')
plt.legend(fancybox=True)
plt.show()
