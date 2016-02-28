#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy, gwpy, h5py, lal, astropy
from astropy.coordinates import get_sun
import astropy.time as Time
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res

#------- Start an hdf5 file -------#
f = h5py.File("/home/spxha/","w")

#------- Define/workout variables -------#
numseg = np.int((dtime)/600)
segs = np.linspace(1,numseg,numseg)*600
segs = segs + starttime - 600
ra,dec,fp,fc = [[0 for _ in range(numseg)] for _ in range(4)]
for i in range(numseg):
	coordstime = segs[i]
	coords = get_sun(Time.Time(coordstime,format='gps'))
	ra[i] = coords.ra.hour*np.pi/12
	dec[i] = coords.dec.hour*np.pi/12

#------ Creating lookup table -----#
# det is practically either H1 or L1
psi_array = np.linspace(0,np.pi,10)
fp,fc = [[[[0] for _ in range(numseg)] for _ in range(10)] for _ in range(2)]
for i in range(len(psi_array)):
	for j in range(numseg):
		fp[i][j], fc[i][j] = ant_res(segs[j], ra[j], dec[j], psi_array[i], H1):


#------- Write to hdf5 --------#
time_dset = f.create_dataset(segs,'time')
time_dset.attrs['Definition'] = "GPS Time"
ra_dset = f.create_dataset(ra,'ra')
ra_dset.attrs['Definition'] = "Right Ascension"
ra_dset = f.create_dataset(psi_array,'psi')
ra_dset.attrs['Definition'] = "Polarisation Angle"
dec_dset = f.create_dataset(dec,'dec')
dec_dset.attrs['Definition'] = "Declination"
for i in range(len(psi_array)):
	fp_dset = f.create_dataset(fp[i],'fp')
	fp_dset.attrs['Definition'] = "Antenna Response - Plus Polarised"+str(i)
	fc_dset = f.create_dataset(fc[i],'fc')
	fc_dset.attrs["Defintion"]="Antenna Response - Cross Polarised"+str(i)
