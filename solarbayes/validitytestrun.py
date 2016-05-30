def validity_test(starttime, endtime, h0_max=0.001):
	#------- Packages ---------#
	import numpy as np
	import astropy, gwpy, h5py, lal, sys, os
	from astropy.coordinates import get_sun
	import astropy.time as Time
	from matplotlib.backends.backend_pdf import PdfPages
	from scipy.signal import butter, filtfilt
	import matplotlib.pyplot as plt
	from gwpy.timeseries import TimeSeries
	from antres import antenna_response as ant_res
	from scipy.misc import logsumexp
	from notchfilt import get_filter_coefs, filter_data
	from optparse import OptionParser
	from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer, AdaptiveETA, AbsoluteETA, AdaptiveTransferSpeed

	widgets = ['Finding Likelihood ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
           ' ', AbsoluteETA(), ' ', AdaptiveTransferSpeed()]

	#-------- Importing, filtering and timeshifting data ----------#
	print 'Reading data'
	durationH = endtime - starttime 	# same for durationL, but not used anyway
	oldstarttime = starttime  			# for use later
	oldendtime = endtime 				# for use later
	Xspacing = 2.44140625E-4
	gpsTime = np.linspace(starttime,endtime,int(1/Xspacing))
	pathtoinput = "/home/spxha/"
	strainH = TimeSeries.read(pathtoinput+'S6framesH1.lcf',channel='H1:LDAS-STRAIN', start=starttime, end=endtime)
	strainL = TimeSeries.read(pathtoinput+'S6framesL1.lcf',channel='L1:LDAS-STRAIN', start=starttime, end=endtime)
	num_points = int(durationH/Xspacing)
	h0_min=0.000001
	h0_vals_num=30
	wm = 'fake_signal'

	# Adding in a fake signal
	# frequency in the middle of the range we're interested in 125 Hz
	omega = 250.0*np.pi
	amplitude = 0.5e-26 # (hopefully) middle of h0 values
	# the signal to introduce is of the form amplitude*np.sin(omega*t)
	# first get timeL and timeH in synch with the sun's position

	print 'Getting the detectors in sync'
	timeH = np.arange(starttime, endtime, Xspacing)
	timeL = timeH
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
		if ((timeL[int(i/Xspacing)]>seg30[i])&(timeL[int(i/Xspacing)]<seg30[i+1])):
			coordstime=seg30[i]
			coords = get_sun(Time.Time(coordstime,format='gps'))
			tdelay[i] = lal.ArrivalTimeDiff(detH1.location, detL1.location, coords.ra.hour*np.pi/12, coords.dec.hour*np.pi/12, tgps)
		else:
			pass
	tdelay = np.repeat(tdelay,int(30/Xspacing))

	# make sure tdelay and timeL are of same length in case integer-ing caused slight inconsistency.
	b = np.ones(len(timeL)-len(tdelay))*tdelay[-1]
	tdelay = np.append(tdelay,b)
	newtimeL = newtimeL - tdelay

	# add in fake signal
	print 'Insert fake signal'
	fakeL = amplitude*np.sin(omega*newtimeL)
	fakeH = amplitude*np.sin(omega*newtimeH)
	newstrainL0=newstrainL[tdelayidx:len(newstrainL)]
	newstrainH0=newstrainH[0:len(newstrainL0)]
	strainH = strainH + fakeH
	strainL = strainL + fakeL

	# H1 and L1 are now in sync and filtered between 100 and 150 Hz.

	#----- Applying a bandpass filter -----#
	print 'Filtering data'
	coefsL = get_filter_coefs('L1')
	coefsH = get_filter_coefs('H1')
	strainL = filter_data(strainL,coefsL)
	strainH = filter_data(strainH,coefsH)

	#----------- Down-sampling ------------#
	print 'Down-sampling'
	Xspacing = Xspacing*32
	num_points = int(durationH/Xspacing)
	newtimeL, newtimeH, newstrainH, newstrainL = [[0 for _ in range(num_points)] for _ in range(4)]
	for i in range(num_points):
		j = 32*i + 16
		newstrainH[i] = np.mean(strainH[j-16:j+16])
		newstrainL[i] = np.mean(strainL[j-16:j+16])
		newtimeH[i] = timeH[j]
		newtimeL[i] = timeL[j]

	newstrainH = newstrainH[76800:len(newstrainH)]
	newstrainL = newstrainL[76800:len(newstrainL)]
	newtimeH = newtimeH[76800:-1]
	newtimeL = newtimeL[76800:-1]
	newtimel = newtimeL
	starttime = starttime + 150
	durationH = newtimeH[0] - newtimeH[len(newtimeH)]
	num_points = int(durationH/Xspacing)
	tdelayidx = int(tdelay / Xspacing)


	############################################################
	#------------ Finding probability distribution ------------#
	#   ------------ Defining some stuff for p -------------   #
	print 'Finding likelihood Part 1/2'
	numseg = int((durationH)/600)
	segs = np.linspace(0,numseg,numseg+1)*600
	segs = segs + newtimeL[0]
	ra,dec,fp,fc = [[0 for _ in range(numseg+1)] for _ in range(4)]
	for i in range(numseg+1):
		coordstime = segs[i]
		coords = get_sun(Time.Time(coordstime,format='gps'))
		ra[i] = coords.ra.hour*np.pi/12
		dec[i] = coords.dec.hour*np.pi/12
	psi_array = np.linspace(0,np.pi,10)
	dpsi = psi_array[1]-psi_array[0]
	sigmaA = 10.0
	h0min = h0_min*np.std(newstrainH)
	h0max = h0_max*np.std(newstrainH)
	h0_array = np.linspace(h0min,h0max,h0_vals_num)
	invSigma0 = np.array([[(1./sigmaA**2), 0.], [0., (1./sigmaA**2)]])
	detSigma0 = sigmaA**4
	dX = newstrainH0
	dY = newstrainL0
	FcX0, FpX0, FcY0, FpY0 = [[0 for _ in range(num_points)] for _ in range(4)]
	for i in range(num_points):
		FpX0[i], FcX0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'H1')
		FpY0[i], FcY0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'L1')
	p = [[0  for _ in range(len(h0_array))] for _ in range(2)]
	ppsi = [0 for _ in range(len(psi_array))]
	logdpsi_2 = np.log(0.5*dpsi)

	cos2pi, sin2pi = [[0 for _ in range(len(psi_array))] for _ in range(2)]
	FpX, FcX, FpY, FcY = [[[0 for _ in range(num_points)] for _ in range(len(psi_array))] for _ in range(4)]
	pbar = ProgressBar(widgets=widgets, max_value=len(psi_array)-1)
    pbar.start()
	for k in range(len(psi_array)):
		cos2pi[k] = np.cos(2*psi_array[k])
		sin2pi[k] = np.sin(2*psi_array[k])
		for i in range(num_points):
			FpX[k][i] = FpX0[i]*cos2pi[k] + FcX0[i]*sin2pi[k]
			FcX[k][i] = FcX0[i]*cos2pi[k] - FpX0[i]*sin2pi[k]
			FpY[k][i] = FpY0[i]*cos2pi[k] + FcY0[i]*sin2pi[k]
			FcY[k][i] = FcY0[i]*cos2pi[k] - FpY0[i]*sin2pi[k]
        pbar.update(k)
    pbar.finish()
    print
	print 'Finding likelihoot Part 2/2. This will take a while... '
	pbar = ProgressBar(widgets=widgets, max_value=num_points-1)
	pbar.start()
	for i in range(num_points):
		d = np.array([dX[i], dY[i]])
		d.shape = (2,1)
		if (i + int(60/Xspacing)<num_points):
			int1 = i + int(60/Xspacing)
		else:
			int1 = i
		if (i - int(60/Xspacing)>0):
			int0 = i - int(60/Xspacing)
		else:
			int0 = 0
		sigmaX = np.std(newstrainH0[int0:int1])
		sigmaY = np.std(newstrainL0[int0:int1])
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
			p[0][j] += logdpsi_2 + logsumexp([logsumexp(ppsi[:-1]), logsumexp(ppsi[1:])])
		pbar.update(i)
	pbar.finish()
	################################
	#--------- Background ---------#
	################################
	background_intervals = np.linspace(1,10,10)*9000
	ra_background, background_tdelay, dec_background, coords_background = [[[0] for _ in range(10)] for _ in range(4)]
	background_intervals = background_intervals.astype(int)
	for j in range(len(background_intervals)):
		coords_background[j]=get_sun(Time.Time(gpsStartH-background_intervals[j],format='gps'))
		ra_background[j]  = coords_background[j].ra.hour  * np.pi/12
		dec_background[j] = coords_background[j].dec.hour * np.pi/12
		background_tdelay[j] = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra_background[j], dec_background[j], tgps)


	tdelayidx = int(background_tdelay[5] / Xspacing)
	newstrainL=newstrainL[tdelayidx:len(newstrainL)]
	if len(newstrainL)<len(newstrainH):
		newstrainH=newstrainH[0:len(newstrainL)]
	numseg = int((durationH)/600)
	segs = np.linspace(0,numseg,numseg+1)*600
	segs = segs + newtimeL[0]
	ra,dec,fp,fc = [[0 for _ in range(numseg+1)] for _ in range(4)]

	for i in range(numseg+1):
		coordstime = segs[i]
		coords = get_sun(Time.Time(coordstime,format='gps'))
		ra[i] = coords.ra.hour*np.pi/12
		dec[i] = coords.dec.hour*np.pi/12
	psi_array = np.linspace(0,np.pi,10)
	dpsi = psi_array[1]-psi_array[0]
	sigmaA = 10.0
	h0min = h0_min*np.std(newstrainH)
	h0max = h0_max*np.std(newstrainH)
	h0_array = np.linspace(h0min,h0max,h0_vals_num)
	invSigma0 = np.array([[(1./sigmaA**2), 0.], [0., (1./sigmaA**2)]])
	detSigma0 = sigmaA**4
	dX = newstrainH
	dY = newstrainL
	FcX0, FpX0, FcY0, FpY0 = [[0 for _ in range(num_points)] for _ in range(4)]
	for i in range(num_points):
		FpX0[i], FcX0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'H1')
		FpY0[i], FcY0[i] = ant_res(gpsTime[int(i*Xspacing/600.)], ra[int(i*Xspacing/600.)], dec[int(i*Xspacing/600.)], 0, 'L1')
	ppsi = [0 for _ in range(len(psi_array))]
	logdpsi_2 = np.log(0.5*dpsi)

	cos2pi, sin2pi = [[0 for _ in range(len(psi_array))] for _ in range(2)]
	FpX, FcX, FpY, FcY = [[[0 for _ in range(num_points)] for _ in range(len(psi_array))] for _ in range(4)]
	for k in range(len(psi_array)):
		cos2pi[k] = np.cos(2*psi_array[k])
		sin2pi[k] = np.sin(2*psi_array[k])
		for i in range(num_points):
			FpX[k][i] = FpX0[i]*cos2pi[k] + FcX0[i]*sin2pi[k]
			FcX[k][i] = FcX0[i]*cos2pi[k] - FpX0[i]*sin2pi[k]
			FpY[k][i] = FpY0[i]*cos2pi[k] + FcY0[i]*sin2pi[k]
			FcY[k][i] = FcY0[i]*cos2pi[k] - FpY0[i]*sin2pi[k]
	pbar = ProgressBar(widgets=widgets, max_value=num_points-1)
	pbar.start()
	for i in range(num_points):
		d = np.array([dX[i], dY[i]])
		d.shape = (2,1)
		if (i + int(60/Xspacing)<num_points):
			int1 = i + int(60/Xspacing)
		else:
			int1 = i
		if (i - int(60/Xspacing)>0):
			int0 = i - int(60/Xspacing)
		else:
			int0 = 0
		sigmaX = np.std(newstrainH[int0:int1])
		sigmaY = np.std(newstrainL[int0:int1])
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
			p[1][j] += logdpsi_2 + logsumexp([logsumexp(ppsi[:-1]), logsumexp(ppsi[1:])])
		pbar.update(i)
	pbar.finish()


	# Write into a file.
	if os.path.exists(wm)==False:
		os.mkdir(wm)
	else:
		pass

	np.savetxt(wm+'/p'+str(oldstarttime)+'.txt',p)
	np.savetxt(wm+'/h0'+str(oldstarttime)+'.txt',h0_array)
