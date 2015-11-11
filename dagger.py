# script to write the dag file
import numpy as np

# first read in
timearray=np.array(np.loadtxt('intersect.txt',dtype='f8'))
StartTimes = timearray[:,0]
EndTimes   = timearray[:,1]

# now to write the .dag file
file = open("mainjob.dag", "w")
for i in range(len(StartTimes)):
	file.write('JOB '+str(i+1)+' job.sub\n')
	file.write('VARS '+str(i+1)+' starttime='+str(StartTimes[i])+', endtime='+str(EndTimes[i])+', outputfile='+str(int(StartTimes[i]))+'\n')
file.close()
