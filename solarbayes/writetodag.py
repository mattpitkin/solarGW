def writetodag(starttime, endtime, h0max='adaptive',timelimit='none'):
	import numpy as np
	import os
	fname1 = str(starttime)
	fname2 = str(endtime)
	if timelimit=='none':
		filepath='../intersect_old.txt'
	elif timelimit=='30 mins':
		filepath='../intersect.txt'
	else:
		print 'Error... illegal timelimit value'
		exit()
	timearray = np.array(np.loadtxt(filepath,dtype='f8'))
	StartTimes = timearray[:,0]
	EndTimes   = timearray[:,1]
	if h0max=='adaptive':
		if 'week' in endtime:
			h0_max = 0.0004
		elif 'month' in endtime:
			h0_max = 0.0001
		elif 'all' in endtime:
			h0_max = 0.00002
		else:
			h0_max = h0max
	elif h0max=='def':
		h0_max = 0.0001
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
		if starttime == 'all':
			starttime = 931076896
			endtime = 971614889

	# Find nearest value to targeted starttime and endtime
	newStartTimes = StartTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]
	newEndTimes   =   EndTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]

	# Now everything is ready, write to a dag file
	file = open("bayesjob"+fname1+'_'+fname2+'_'+timelimit+".dag", "w")
	for i in range(len(newStartTimes)):
		file.write('JOB '+str(i+1)+' bayesjob.sub\n')
		file.write('VARS '+str(i+1)+' starttime="'+str(int(newStartTimes[i]))+'" endtime="'+str(int(newEndTimes[i]))+'" h0="'+str(h0_max)+'" Proc="'+str(i+1)+'"'+'\n')
	file.close()
