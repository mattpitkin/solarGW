# code to add up the probability for a month or a week, prints the duration, and plots the probability distribution
def plotdist(starttime,endtime):
	# starttime and endtime should be in format 'week'+str(i) or 'month'+str(i) only.
	# In the near future it will take 'all' as starttime and endtime to plot the entire data.
	import numpy as np
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.pyplot as plt
	import os
	fname =	'probdist_cumulative.pdf'
	if os.path.exists('/home/spxha/')==True:
		pathtointersect = '/home/spxha/solarGW/intersect_old.txt'
	else:
		pathtointersect = '../../../intersect_old.txt'
	timearray = np.array(np.loadtxt(pathtointersect,dtype='f8'))
	StartTimes = timearray[:,0]
	EndTimes   = timearray[:,1]
	if starttime=='all' and endtime=='all':
		starttime = 931076896.0
		endtime   = 971614889.0
	duration = 0
	newStartTimes = StartTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes-endtime)).argmin()]
	newEndTimes   =   EndTimes[(np.abs(StartTimes-starttime)).argmin():(np.abs(EndTimes-endtime)).argmin()]
	p_array,h0_array = [[[0 for _ in range(30)] for _ in range(len(newStartTimes))] for _ in range(2)]
	j = 0
	for i in range(len(newStartTimes)):
		pathi = 'p'+str(int(newStartTimes[i]))+'.txt'
		if os.path.exists(pathi)==True:
			j +=1
			print i
			p_array[i]  = np.array(np.loadtxt('p'+str(int(newStartTimes[i]))+'.txt',dtype='float'))
			h0_array[i] = np.array(np.loadtxt('h0'+str(int(newStartTimes[i]))+'.txt',dtype='float'))
			duration += newEndTimes[i]-newStartTimes[i] - 150
		else:
			pass
	print j
	p_sum_array = [0.0 for _ in range(30)]
	p_sum_array = np.array(p_sum_array)
	h0_array = np.array(h0_array)
	print p_array[1]
	p_array = np.array(p_array)
	for i in range(len(p_array)-1):
		p_sum_array += p_array[i+1]
	h0_mean_array = h0_array.mean(0)
	p = np.exp(p_sum_array-np.max(p_sum_array))
	print h0_mean_array,p

	# Making a cumulative plot to get the 95% value.
	prob = [0 for _ in range(len(p))]
	for i in range(30):
		d = [0 for _ in range(i+1)]
		print i
		for j in range(i+1):
			print i,j
			d[j] = p[i-j]
			print d[j]
		prob[i] = np.sum(d)
		print prob[i]
	prob = prob/(np.max(prob))
	print prob
	# plot the distribution
	with PdfPages(fname) as pdf:
		fig1 = plt.figure()
		plt.plot(h0_mean_array,prob,'+')
#		plt.title('Cumulative probability for '+str(int(duration/3600))+' hours')
		plt.xlabel(r'$h_0$')
		plt.ylabel('Normalised cumulative probability')
		pdf.savefig(fig1)
		plt.close()
	print (h0_mean_array[10]+h0_mean_array[11])/2
