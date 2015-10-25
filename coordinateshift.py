def get_coordsdiff(T1=967966720):
	from astropy.coordinates import get_sun
	import astropy.time as Time
	
	T1 =  967966720
	t = Time.Time(T1,format='gps')
	t.coords = get_sun(t)
	coords = get_sun(t)
	T2 = T1 + 2400
	t2 = Time.Time(T2,format='gps')
	coords2 = get_sun(t2)
	dec = coords.dec.hour
	dec2 = coords2.dec.hour
	
	# Difference of declination after 40 minutes (duration of 1 hdf5 file)
	diff_dec = dec2-dec
	print('Difference in declination between T = 967966720 GPS time and 40 minutes later',
	'is', diff_dec)
	# Doing the same for right ascension:
	diff_ra = coords2.ra.hour - coords.ra.hour
	print('Difference in right ascention between T = 967966720 GPS time and 40 minutes later',
	'is', diff_ra)