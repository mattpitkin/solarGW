#!/usr/bin/env python

# plotting all the data
from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
import numpy as np
import gwpy
import glue


data = np.loadtxt('/home/spxha/solarGW/intersect6.txt',dtype='f8')
TimesStart = data[:,0]
TimesEnd = data[:,1]
Xspacing = 2.44140625E-4
for i in range(len(TimesStart)):
	print i	
	duration = TimesEnd[i]-TimesStart[i]
	timeseriesi = np.linspace(TimesStart[i], TimesEnd[i], int(duration/Xspacing))
	straini = TimeSeries.read('/home/spxha/S6framesH1.lcf',channel='H1:PSL-ODC_CHANNEL_OUT_DQ',start=TimesStart[i], end=TimesEnd[i], format='lcf')

	if i==0:
		strain = straini
		timeseries = timeseriesi
	else:
		strain = np.append(strain, straini)
		timeseries = np.append(timeseries, timeseriesi)
	del duration, timeseriesi, straini

plt.figure()
plt.plot(timeseries,strain)
plt.savefig('/home/spxha/Timeseries.png')
