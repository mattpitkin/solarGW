# plotting all the data
from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('intersect6.txt',dtype='f8')
TimesStart = data[:,0]
TimesEnd = data[:,1]
Xspacing = 2.44140625E-4
for i in range(len(TimesTimeseries)):
	duration[i] = TimesEnd[i]-TimesStart[i]
	timeseriesi = np.linspace(TimesStart[i], TimesEnd[i], int(duration[i]/Xspacing))
	straini = TimeSeries.read('S6framesH1.lcf',channel='H1:LDAS-STRAIN',start=TimesStart[i], end=TimesEnd)

	if i==0:
		strain = straini
		timeseries = timeseriesi
	else:
		strain = np.append(strain, straini)
		timeseries = np.append(timeseries, timeseriesi)

plt.figure()
plt.plot(timeseries,strain)
plt.savefig('/home/spxha/TimeSeries.png')
