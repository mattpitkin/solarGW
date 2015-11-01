from astropy.coordinates import get_sun
import astropy.time as Time
import numpy as np
import lal
import matplotlib.pyplot as plt

T0 =  967966720
detMap = {'H1': lal.LALDetectorIndexLHODIFF, 'L1':
lal.LALDetectorIndexLLODIFF}
detH1 = lal.CachedDetectors[detMap['H1']]
detL1 = lal.CachedDetectors[detMap['L1']]
t,T1,coords,dec,ra,tdelay,tgps = [[[0] for _ in range(100000)] for _ in range(7)]
for i in range(100000):
	T1[i] = T0 + i
	t[i] = Time.Time(T1[i],format='gps')
	t[i].coords = get_sun(t[i])
	coords[i] = get_sun(t[i])
	dec[i] = coords[i].dec.hour * np.pi/12
	ra[i] = coords[i].ra.hour *np.pi/12
	tgps[i] = lal.LIGOTimeGPS(T1[i], 0)
	tdelay[i] = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra[i], dec[i], tgps[i])
plt.plot(T1, tdelay)
plt.show()
