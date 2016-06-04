# simple script to make sure there are no segments larger than 30 minutes.

import numpy as np
timearray = np.array(np.loadtxt('intersect.txt',dtype='f8'))
StartTimes = timearray[:,0]
EndTimes   = timearray[:,1]
for i in range(len(StartTimes)):
	if EndTimes[i]-StartTimes[i]>1800:
		for j in range(int((EndTimes[i]-StartTimes[i])/1800)+1):
			StartTimes = np.append(StartTimes,StartTimes[i]+1800*j)
			EndTimes = np.append(EndTimes, StartTimes[i]+1800*j-1)
	if EndTimes[i]-StartTimes[i]<250:
		StartTimes = np.delete(StartTimes,[StartTimes[i]])
		EndTimes   = np.delete(EndTimes,[EndTimes[i]])
newStartTimes = np.sort(StartTimes)
newEndTimes   = np.sort(EndTimes)
newEndTimes   = newEndTimes.astype(int)
newStartTimes = newStartTimes.astype(int)
np.savetxt('newintersect.txt',np.c_[newStartTimes,newEndTimes], fmt='%i')
