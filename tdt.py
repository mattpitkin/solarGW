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
t,T1,coords,dec,ra,tdelay,tgps = [[[0] for _ in range(100000)] for _ in range(7)] # initialising arrays
T1_pre = np.linspace(1,100000,100000)
# looping here because lal.LIGOTimeGPS does not seem to like numpy arrays.
for i in range(100000):
	T1[i] = T0 + i
# T1 = T0 + T1_pre
t = Time.Time(T1,format='gps')
coords = get_sun(t)
dec = coords.dec.hour * np.pi/12
ra  = coords.ra.hour  * np.pi/12
# again looping for same reason as above
# for i in range(100000):
# 	tgps[i] = lal.LIGOTimeGPS(T1[i], 0)
# 	tdelay[i] = lal.ArrivalTimeDiff(detH1.location, detL1.location, ra[i], dec[i], tgps[i])
# j_s= np.linspace(1,162,162)*600
# ddt = [[0] for _ in range(162)]
# a = np.linspace(1,160,160)*600
# a = a.astype(int)
# for j in a:
# 	ddt[j/600]=tdelay[j+600]-tdelay[j]
# T_s = T0 + j_s

j_s2= np.linspace(1,320,320)*300
ddt2 = [[0] for _ in range(320)]
a2 = np.linspace(1,315,315)*300
a2 = a2.astype(int)
for j in a2:
	ddt2[j/300]=tdelay[j+300]-tdelay[j]
T_s2 = T0 + j_s2

j_s3= np.linspace(1,1600,1600)*60
ddt3 = [[0] for _ in range(1600)]
a3 = np.linspace(1,1500,1500)*60
a3 = a3.astype(int)
for j in a3:
	ddt3[j/60]=tdelay[j+60]-tdelay[j]
T_s3 = T0 + j_s3

j_s4= np.linspace(1,3200,3200)*30
ddt4 = [[0] for _ in range(3200)]
a4 = np.linspace(1,3100,3100)*30
a4 = a4.astype(int)
for j in a4:
	ddt4[j/30]=tdelay[j+30]-tdelay[j]
T_s4 = T0 + j_s4

# plt.plot(T_s[1:160], ddt[1:160],color="blue",label='10 minute segments')
plt.plot(T_s2[1:315], ddt2[1:315],color="red",label='5 minute segments')
plt.plot(T_s3[1:1500], ddt3[1:1500],color="green",label='1 minute segments')
plt.plot(T_s4[1:3100], ddt4[1:3100],color="black",label='1/2 minute segments')
plt.xlabel('GPS Time (s)')
plt.ylabel('Change in time delay')
# plt.title('Change in time delay for different time segments')
plt.legend(fancybox=True)
plt.show()
