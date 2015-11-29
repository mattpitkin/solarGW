import numpy as np
import os, h5py
from upperlim import uplim

# Read in the intersect6.txt file
timearray  = np.array(np.loadtxt('intersect6.txt',dtype='f8'))
StartTimes = timearray[:,0]
EndTimes   = timearray[:,1]
duration = 0
strain = 0
for i in range(len(StartTimes)):
	h5path =  'scaratch/spxha/'+str(StartTimes[i])+'_'+str(EndTimes[i])
	if os.path.exists(h5path):
		durationi = EndTimes[i] - StartTimes[i]
		duration = duration + durationi
		f = h5py.File(h5path+str(StartTimes[i])+'.hdf5','r')
		for j in range(len(f['Signal'])):
			strainj = f['Signal/Signal'+str(j)].value
			if j==0:
				strain = strain + strainj
result = uplim(strain,duration)
print result
