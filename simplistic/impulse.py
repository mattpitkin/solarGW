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
from scipy.signal import butter
from scipy.signal import filtfilt

timeseries = np.linspace(0,600,600)
strain = np.zeros(600)
strain[0]=1
ord = 4
Wn = [100.0/2918.0,300.0/2918.0]
type = 'bandpass'
b,a = butter(ord,Wn, btype=type)
strain = filtfilt(b,a,strain)
plt.plot(timeseries,strain)
plt.xlabel('Time (s)')
plt.ylabel('Strain')
plt.title('Strain data from a filtered impulse')
plt.show()
