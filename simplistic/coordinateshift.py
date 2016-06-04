def get_coordsdiff(x,T1=967966720):
	from astropy.coordinates import get_sun
	import astropy.time as Time
	import numpy as np
	import lal


	T1 =  967966720
	t = Time.Time(T1,format='gps')
	t.coords = get_sun(t)
	coords = get_sun(t)
	T2 = T1 + x
	t2 = Time.Time(T2,format='gps')
	coords2 = get_sun(t2)
	dec = coords.dec.hour * np.pi/12
	ra = coords.ra.hour *np.pi/12
	ra2 = coords.ra.hour * np.pi/12
	dec2 = coords2.dec.hour * np.pi/12

	detMap = {'H1': lal.LALDetectorIndexLHODIFF, 'L1':
	lal.LALDetectorIndexLLODIFF}
	detH1 = lal.CachedDetectors[detMap['H1']]
	detL1 = lal.CachedDetectors[detMap['L1']]
	tgps = lal.LIGOTimeGPS(T1, 0)
	tgps2 = lal.LIGOTimeGPS(T2, 0)

	coords = get_sun(Time.Time(967966720.0,format='gps'))
	ra  = coords.ra.hour  * np.pi/12
	dec = coords.dec.hour * np.pi/12


	tdelay1 = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra, dec, tgps)
	tdelay2 = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra2, dec2, tgps2)
	diff_dec = tdelay2-tdelay1

