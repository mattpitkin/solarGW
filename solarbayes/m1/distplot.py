# code to add up the probability for a month or a week, prints the duration, and plots the probability distribution
def plotdist(starttime,endtime):
	# starttime and endtime should be in format 'week'+str(i) or 'month'+str(i) only.
	# In the near future it will take 'all' as starttime and endtime to plot the entire data.
	import numpy as np
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.pyplot as plt
	import os
	fname1 = str(starttime)
	fname2 = str(endtime)
	timearray = np.array(np.loadtxt('../intersect.txt',dtype='f8')
	StartTimes = timearray[:,0]
	EndTimes   = timearray[:,1]
	for i in range(70):
		if starttime == 'month'+str(i):
			starttime = 931076896 + 2592000 * i
			wm = 'm'
		else:
			pass
		if endtime   == 'month'+str(i):
			endtime = 931076896 + 2592000 *i
			wm = 'm'
		else:
			pass
		if starttime == 'week'+str(i):
			starttime = 931076896 + 604800 * i
			wm = 'w'
		else:
			pass
		if endtime   == 'week'+str(i):
			endtime   = 931076896 + 604800 * i
			wm = 'w'
		else:
			pass
	newStartTimes = StartTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]
	newEndTimes   =   EndTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes  -  endtime)).argmin()]
	p_array,h0_array = [[[0 for _ in range(30)] for _ in range(len(newStartTimes))] for _ in range(2)]
	for i in range(len(newStartTimes)):
		p_array[i]  = np.array(np.loadtxt('p'+str(int(newStartTimes[i]))+wm,dtype='float'))
		h0_array[i] = np.array(np.loadtxt('h0'+str(int(newStartTimes[i]))+wm,dtype='float'))

	for i in range(len(p_array)-1):
		p_sum_array += p_array[i+1]
		h0_mean_array = h0_array.mean(0)

	# plot the distribution
	with PdfPages(fname) as pdf:
		fig1 = plt.figure()
		plt.plot(h0_mean_array,p_sum_array,'+')
		plt.title('Probability Distribution for GPS'+startname+'_'+endtime)
		plt.xlabel('h0 '+str(h0_min)+' '+str(h0_max))
		plt.ylabel('p')
		pdf.savefig(fig1)
		plt.close()
#	print 'Plot saved as', fname, 'Did we find GWs?!'
