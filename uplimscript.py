import numpy as np
import os, h5py
from upperlim import uppperlim

# Read in the intersect6.txt file
timearray  = np.array(np.loadtxt('intersect6.txt',dtype='f8'))
StartTimes = timearray[:,0]
EndTimes   = timearray[:,1]
duration = 1
strain = 0
for i in range(len(StartTimes)-1):
	h5path =  'scratch/spxha/'+str(StartTimes[i])+'_'+str(EndTimes[i])
	if os.path.exists(h5path):
		durationi = EndTimes[i] - StartTimes[i]
		if (durationi>180):
			duration = duration + durationi - 150
			f = h5py.File(h5path+str(StartTimes[i])+'.hdf5','r')
			for j in range(5,len(f['Signal'])-1):
				strainj = f['Signal/Strail'+str(j)].value
				if j==0:
					strain = strain + strainj
				else:
					pass
		else:
			pass
	else:
		pass
result = upperlim(strain,duration)
