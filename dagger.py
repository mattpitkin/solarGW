# script to write the dag file
import numpy as np

# first read in
timearray=np.array(np.loadtxt('intersect6.txt',dtype='f8'))
StartTimes = timearray[:,0]
EndTimes   = timearray[:,1]

# now to write the .dag file
file = open("mainjob.dag", "w")
for i in range(len(StartTimes)):
	file.write('JOB '+str(i+1)+' testjob.sub\n')
	file.write('VARS '+str(i+1)+' starttime="'+str(int(StartTimes[i]))+'" endtime="'+str(int(EndTimes[i]))+'"'+'\n')
file.close()
