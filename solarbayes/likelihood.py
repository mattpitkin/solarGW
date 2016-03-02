def likelihood(starttime=969062862, endtime=969063609, h0_factor=3, h0_vals_num=25):
	#------- Packages ---------#
	import numpy as np
	import astropy, gwpy, h5py, lal
	from astropy.coordinates import get_sun
	import astropy.time as Time
	from scipy.signal import butter
	from scipy.signal import filtfilt
	from gwpy.timeseries import TimeSeries
	from antres import antenna_response as ant_res
	from scipy.misc import logsumexp
	from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer, AdaptiveETA, AbsoluteETA, AdaptiveTransferSpeed

	widgets = ['Finding Likelihood ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
           ' ', AbsoluteETA(), ' ', AdaptiveTransferSpeed()]

	################################################################
	#-------- Importing, filtering and timeshifting data ----------#
	print 'Reading data'
	gpsStartH = starttime
	durationH = endtime - starttime
	gpsEndH   = endtime
	gpsStartL = gpsStartH
	durationL = durationH
	gpsEndL   = gpsEndH
	Xspacing = 2.44140625E-4
	gpsTime = np.linspace(starttime,endtime,int(1/Xspacing))
	pathtoinput = "/home/spxha/"
	strainH = TimeSeries.read(pathtoinput+'S6framesH1.lcf',channel='H1:LDAS-STRAIN', start=starttime, end=endtime)
	strainL = TimeSeries.read(pathtoinput+'S6framesL1.lcf',channel='L1:LDAS-STRAIN', start=starttime, end=endtime)

	#---------------------------
	# Applying a bandpass filter
	#---------------------------
	print 'Filtering data'
	ord = 4
	Wn = [100.0/2918.0,150.0/2918.0]
	type = 'bandpass'
	bL,aL = butter(ord,Wn, btype=type)
	bH,aH = butter(ord,Wn, btype=type)
	strainL = filtfilt(bL,aL,strainL)
	strainH = filtfilt(bH,aH,strainH)

	print 'Getting the detectors in sync'
	timeH = np.arange(gpsStartH, gpsEndH, Xspacing)
	timeL = np.arange(gpsStartL, gpsEndL, Xspacing)
	timel=timeL
	detMap = {'H1': lal.LALDetectorIndexLHODIFF, 'L1':
	lal.LALDetectorIndexLLODIFF}
	detH1 = lal.CachedDetectors[detMap['H1']]
	detL1 = lal.CachedDetectors[detMap['L1']]
	tgps = lal.LIGOTimeGPS(gpsStartH, 0)

	#---------- Get right ascension and declination of source in radians ----------#
	numseg30 = int((endtime-starttime)/30.)
	seg30 = gpsStartH + 30*np.linspace(1,numseg30,numseg30) # 30 second update rate
	tdelay = [[0] for _ in range(numseg30)]
	for i in range(numseg30-1):
		if ((timel[int(i/Xspacing)]>seg30[i])&(timel[int(i/Xspacing)]<seg30[i+1])):
			coordstime=seg30[i]
			coords = get_sun(Time.Time(coordstime,format='gps'))
			tdelay[i] = lal.ArrivalTimeDiff(detH1.location, detL1.location, coords.ra.hour*np.pi/12, coords.dec.hour*np.pi/12, tgps)
		else:
			pass
	tdelay = np.repeat(tdelay,int(30/Xspacing))

	# make sure tdelay and timeL are of same length in case integer-ing caused slight inconsistency.
	b = np.ones(len(timeL)-len(tdelay))*tdelay[-1]
	tdelay = np.append(tdelay,b)

	timeL = timeL - tdelay
	# H1 and L1 are now in sync and filtered between 100 and 150 Hz.

	############################################################
	#------------ Finding probability distribution ------------#
	#------- Defining some stuff for p ------#
	print 'Finding likelihood Part 1/2'
	numseg = int((durationH)/600)
	segs = np.linspace(0,numseg,numseg+1)*600
	segs = segs + starttime
	ra,dec,fp,fc = [[0 for _ in range(numseg+1)] for _ in range(4)]
	for i in range(numseg+1):
		coordstime = segs[i]
		coords = get_sun(Time.Time(coordstime,format='gps'))
		ra[i] = coords.ra.hour*np.pi/12
		dec[i] = coords.dec.hour*np.pi/12
	psi_array = np.linspace(0,np.pi,10)
	dpsi = psi_array[1]-psi_array[0]
	sigmaA = 10.0
	h0_array = np.linspace(0.001*np.std(strainH),h0_factor*np.std(strainH),h0_vals_num)
	invSigma0 = np.array([[(1./sigmaA**2), 0.], [0., (1./sigmaA**2)]])
	detSigma0 = sigmaA**4
	dX = strainH
	dY = strainL
	FcX0, FpX0, FcY0, FpY0 = [[0 for _ in range(int(durationH/Xspacing))] for _ in range(4)]
	for i in range(int(durationH/Xspacing)):
		FpX0[i], FcX0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'H1')
		FpY0[i], FcY0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'L1')
	p = [0  for _ in range(len(h0_array))]
	ppsi = [0 for _ in range(len(psi_array))]
	logdpsi_2 = np.log(0.5*dpsi)

	cos2pi, sin2pi = [[0 for _ in range(psi_array)] for _ in range(2)]
	FpX, FcX, FpY, FcY = [[[ ]]]
	for k in range(len(h0_array)):
		cos2pi[k] = np.cos(2*psi_array[k])
		sin2pi[k] = np.sin(2*psi_array[k])
		for i in range():
			FpX[k][i] = FpX0[i]*cos2pi[k] + FcX0[i]*sin2pi[k]
			FcX[k][i] = FcX0[i]*cos2pi[k] - FpX0[i]*sin2pi[k]
			FpY[k][i] = FpY0[i]*cos2pi[k] + FcY0[i]*sin2pi[k]
			FcY[k][i] = FcY0[i]*cos2pi[k] - FpY0[i]*sin2pi[k]

	print 'Finding likelihoot Part 2/2. This will take a while... '
	pbar = ProgressBar(widgets=widgets, max_value=int(durationH/Xspacing)-1)
	pbar.start()
	for i in range(int(durationH/Xspacing)):
		d = np.array([dX[i], dY[i]])
		d.shape = (2,1)
		if (i + int(60/Xspacing)<int(durationH/Xspacing)):
			int1 = i + int(60/Xspacing)
		else:
			int1 = i
		if (i - int(60/Xspacing)>0):
			int0 = i - int(60/Xspacing)
		else:
			int0 = 0
		sigmaX = np.std(strainH[int0:int1])
		sigmaY = np.std(strainL[int0:int1])
		C = np.array([[sigmaX**2, 0.], [0., sigmaY**2]])
		invC = np.array([[(1./sigmaX**2), 0.], [0., (1/sigmaY**2)]])
		detC = sigmaX**2 * sigmaY**2
		for j in range(len(h0_array)):
			for k in range(len(psi_array)):
				M = h0_array[j]*np.array([[FpX[k][i], FpY[k][i]], [FcX[k][i], FcY[k][i]]])
				M = np.array([[M[0][0][0],M[0][1][0]],[M[1][0][0], M[1][1][0]]])
				invSigma = np.dot(M.T, np.dot(invC, M)) + invSigma0
				Sigma = np.linalg.inv(invSigma)
				detSigma = np.linalg.det(Sigma)
				chi = np.dot(Sigma, np.dot(M.T, np.dot(invC, d)))
				ppsi[k]    = 0.5*np.log(detSigma) - 0.5*np.log(16.*np.pi**4*detSigma0*detC) -  0.5*(np.vdot(d.T, np.dot(invC, d)) + np.vdot(chi.T, np.dot(invSigma, chi)))
			p[j] += logdpsi_2 + logsumexp([logsumexp(ppsi[:-1]), logsumexp(ppsi[1:])])
		pbar.update(i)
	pbar.finish()
	print

	# Write into a file
	f = h5py.File('/home/spxha/probability.hdf5','w')
	f.create_dataset(h0_array,'h0')
	f.create_dataset(p,'p')

	#------ plot the probability distribution
	print 'Producing Plot'
	fname = 'probdist.pdf'
	with PdfPages(fname) as pdf:
		fig1 = plt.figure()
		plt.plot(h0_array,p,'+')
		plt.title('Probability Distribution')
		pdf.savefig(fig1)
		plt.close()
	print 'Plot saved as', fname, 'Did we find GWs?!'
