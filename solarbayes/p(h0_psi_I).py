#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy
import gwpy
import h5py
import lal
from astropy.coordinates import get_sun
import astropy.time as Time
import astropy
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res
from scipy.integrate import quad, dblquad

# we already know fp, fc and psi(?,integration) by now, and obviously h(L) and h(Y)

# need to define standard deviations: sig_x, sig_y, sig_A

# for i in range(time/Xspacing):

V =  h0**2 * sig_Y * ((fc_X/ sig_X**2) + (fc_Y/sig_Y**2))
N1 = - sig_A**2 * (2 * d_X * d_Y * (fp_X * p_Y + fc_X * fc_Y) * h0**2 + d_Y * d_X**2 *h0 * (fp_X**2 + fc_X**2) + sig_A**2 * sig_X**2) )
N2 = 2*(-2*h0**4 * fc_X * fc_Y * fp_X * fp_Y + ((h0**2 fc_X**2*(fp_X**2*h0**2+sig_X**2*sig_A**2)) + (h0**2 fc_Y**2*(fp_Y**2*h0**2+sig_Y**2*sig_A**2))) + sig_A**2 * ((fc_X**2*h0**2*sig_X**2) + (fc_Y**2*h0**2*sig_Y**2) + sig_A**2 *sig_X**2))
p = (1/ (2 * np.pi * sig_A**2 * sig_x[i] * np.sqrt(sig_A**2 + V) ))* np.exp(N1/N2)
d_X**2 *h0 * (fp_X**2 + fc_X**2) + sig_A**2 * sig_X**2)

#def p(Ac,Ap):
#	return  (1/2*np.pi*sig_x*sig_y)*exp()


#ans, err = dblquad(integrand,-np.inf,+np.inf,
#lambda x: -np.inf,
#lambda x: +np.inf,
#args=(,,,))
