import numpy as np
import os, h5py
import matplotlib.pyplot as plt
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer, AdaptiveETA, AbsoluteETA, AdaptiveTransferSpeed
from matplotlib.backends.backend_pdf import PdfPages
from upperlim import upperlim

widgets = ['Reading HDF5s ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
           ' ', AbsoluteETA(), ' ', AdaptiveTransferSpeed()]

# Read in the intersect6.txt file, initialise arrays
timearray  = np.array(np.loadtxt('intersect.txt',dtype='int'))
StartTimes = timearray[:,0]
EndTimes   = timearray[:,1]
duration = 1
strain = 0
backstrain, backj = [[[0] for _ in range(100)] for _ in range(2)]
# Get upper limit
pbar = ProgressBar(widgets=widgets, maxval=4000)
pbar.start()
for i in range(3173):
	h5path =  '/scratch/spxha/'+str(StartTimes[i])+'_'+str(EndTimes[i])
	fullpath=h5path+'/'+str(StartTimes[i])+'.hdf5'
	if os.path.exists(fullpath)==True:
		durationi = EndTimes[i] - StartTimes[i]
		if (durationi>180):
			duration = duration + durationi - 150
			f = h5py.File(fullpath,'r')
			for j in range(5,len(f['Signal'])-1):
				strainj = f['Signal/Strain'+str(j)].value
				strain = strain + strainj
			for i in range(100):
				for j in range(5,len(f['Signal'])-1):
					backj[i] = f['Background'+str(i)+'/Strain'+str(j)].value
					backstrain[i] = backstrain[i] + backj[i]
		else:
			pass
	else:
		pass
	pbar.update(i)
pbar.finish()
print
theupperlimitthativebeenlookingfor = upperlim(strain,duration)
print 'The Upper Limit is ...', theupperlimitthativebeenlookingfor
print 'Duration = ', duration/86400.0, 'days'
print 'Now creating final plot...'
# plot the foreground/background
fname = 'strain.pdf'
with PdfPages(fname) as pdf:
	fig1 = plt.figure()
	# plot the signal as a vertical line
	plt.axvline(strain) #done
	# plot the background as a histogram
	plt.hist(backstrain) #check?
	plt.title('Comparing Strain Signal from the Sun with the Background')
	pdf.savefig(fig1)
	plt.close()
print 'Plot saved as', fname, 'Did we find anything?!'
