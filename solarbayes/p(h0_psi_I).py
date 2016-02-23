#------- Packages ---------#
from __future__ import division
import numpy as np
import astropy, gwpy, h5py, lal, astropy
from astropy.coordinates import get_sun
import astropy.time as Time
from gwpy.timeseries import TimeSeries
from antres import antenna_response as ant_res

# from scipy.integrate import quad, dblquad

# need to define standard deviations: sig_x, sig_y, sig_A

Xspacing = 2.44140625E-4
# for i in range(time/Xspacing):

# data vector (for detectors X and Y)
d = np.array([dX, dY])

# design matrix
M = h0*np.array([[FpX, FpY], [FcX, FcY]])

# the inverse data covariance matrix and the determinant of the covariance matrix
C = np.array([[sigmaX**2, 0.], [0., sigmaY**2]])
invC = np.array([[(1./sigmaX**2), 0.], [0., (1/sigmaY**2)]])
detC = sigmaX**2 * sigmaY**2

# inverse of the amplitude prior covariance matrix
invSigma0 = np.array([[(1./sigmaA**2) 0.], [0., (1./sigmaA**2)]])
detSigma0 = sigmaA**4

# inverse of Sigma
invSigma = np.dot(M.T, np.dot(invC, M)) + invSigma0
Sigma = np.linalg.inv(invSigma)
detSigma = np.linalg.det(Sigma)

# chi
chi = np.dot(Sigma, np.dot(M.T, np.dot(invC, d)))

# log likelihood
p = 0.5*np.log(detSigma) - 0.5*log(16.*np.pi**4*detSigma0*detC) -
    0.5*(np.vdot(d.T, np.dot(invC, d)) + np.vdot(chi.T, np.dot(invSigma, chi)))

