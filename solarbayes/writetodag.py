def writetodag(starttime, endtime, h0max='adaptive'):
	import numpy as np
	timearray = np.array(np.loadtxt('intersect6.txt',dtype='f8'))
	StartTimes = timearray[:,0]
	EndTimes   = timearray[:,1]
	for i in range(70):
		if starttime = 'month'+str(i):
			starttime = 931076896 + 2592000 * i
		else:
			pass
		if endtime = 'month'+str(i):
			endtime = 931076896 + 2592000 *i
		else:
			pass
		if starttime = 'week'+str(i):
			starttime = 931076896 + 604800 * i
		else:
			pass
		if endtime='week'+str(i):
			endtime   = 931076896 + 604800 * i
		else:
			pass

	# Find nearest value to targeted starttime and endtime
	newStartTimes = StartTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]
	newEndTimes   = StartTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]

	# Now everything is ready, write to a dag file
	file = open("bayesjob.dag", "w")
	for i in range(len(newStartTimes)):
		file.write('JOB '+str(i+1)+' bayesjob.sub\n')
		file.write('VARS '+str(i+1)+' starttime="'+str(int(newStartTimes[i]))+'" endtime="'+str(int(newEndTimes[i]))+'"'+'\n')
	file.close()
