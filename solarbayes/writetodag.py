def writetodag(starttime, endtime, h0max='adaptive'):
	import numpy as np
	import os
	timearray = np.array(np.loadtxt('../intersect.txt',dtype='f8'))
	StartTimes = timearray[:,0]
	EndTimes   = timearray[:,1]
	if h0max=='adaptive':
		if 'week' in endtime:
			h0_max = 0.0001
		elif 'month' in endtime:
			h0_max = 0.00001
		else:
			h0_max = h0max
	for i in range(70):
		if starttime == 'month'+str(i):
			starttime = 931076896 + 2592000 * i
		else:
			pass
		if endtime   == 'month'+str(i):
			endtime = 931076896 + 2592000 *i
		else:
			pass
		if starttime == 'week'+str(i):
			starttime = 931076896 + 604800 * i
		else:
			pass
		if endtime   == 'week'+str(i):
			endtime   = 931076896 + 604800 * i
		else:
			pass


	# Find nearest value to targeted starttime and endtime
	newStartTimes = StartTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]
	newEndTimes   =   EndTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]

	# Now everything is ready, write to a dag file
	file = open("bayesjob.dag", "w")
	for i in range(len(newStartTimes)):
		file.write('JOB '+str(i+1)+' bayesjob.sub\n')
		file.write('VARS '+str(i+1)+' starttime="'+str(int(newStartTimes[i]))+'" endtime="'+str(int(newEndTimes[i]))+'" h0="'+str(h0_max)+'"'+'\n')
	file.close()
