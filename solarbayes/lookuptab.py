#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy
import gwpy
import h5py
import lal
from astropy.coordinates import get_sun
import astropy.time as Time
import astropy
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res

#------- Start an hdf5 file -------#
f = h5py.File("/home/spxha/","w")

#------- Define/workout variables -------#
#for every 10 minutes, I was to create (sub-?)arrays of ra, dec,
starttime=931219808

#endime =
numseg = np.int((endtime-starttime)/600)
segs = np.linspace(1,numseg,numseg)*600
segs = segs + starttime - 600
ra,dec,fp,fc = [[0 for _ in range(numseg)] for _ in range(4)]
for i in range(numseg):
	coordstime = segs[i]
	coords = get_sun(Time.Time(coordstime,format='gps'))
	ra[i] = coords.ra.hour*np.pi/12
	dec[i] = coords.dec.hour*np.pi/12

#------ Creating lookup table -----#
# psi = Jeffrey's prior?
# det is practically either H1 or L1
fp, fc = ant_res(gpsTime, ra, dec, psi, det):
