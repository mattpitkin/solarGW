def antenna_response( gpsTime, ra, dec, psi, det ):
  from types import StringType, FloatType
  import lal
  import numpy as np
  from types import StringType, FloatType

  # create detector-name map
  detMap = {'H1': lal.LALDetectorIndexLHODIFF, \
            'H2': lal.LALDetectorIndexLHODIFF, \
            'L1': lal.LALDetectorIndexLLODIFF, \
            'G1': lal.LALDetectorIndexGEO600DIFF, \
            'V1': lal.LALDetectorIndexVIRGODIFF, \
            'T1': lal.LALDetectorIndexTAMA300DIFF, \
            'AL1': lal.LALDetectorIndexLLODIFF, \
            'AH1': lal.LALDetectorIndexLHODIFF, \
            'AV1': lal.LALDetectorIndexVIRGODIFF}

  try:
    detector=detMap[det]
  except KeyError:
    raise ValueError, "ERROR. Key %s is not a valid detector name." % (det)

  # get detector
  detval = lal.CachedDetectors[detector]
  response = detval.response

  # check if gpsTime is just a float or int, and if so convert into an array
  if isinstance(gpsTime, float) or isinstance(gpsTime, int):
    gpsTime = np.array([gpsTime])
  else: # make sure it's a numpy array
    gpsTime = np.copy(gpsTime)

  # if gpsTime is a list of regularly spaced values then use ComputeDetAMResponseSeries
  if len(gpsTime) == 1 or np.unique(np.diff(gpsTime)).size == 1:
    gpsStart = lal.LIGOTimeGPS( gpsTime[0] )
    dt = 0.
    if len(gpsTime) > 1:
      dt = gpsTime[1]-gpsTime[0]
    fp, fc = lal.ComputeDetAMResponseSeries(response, ra, dec, psi, gpsStart, dt, len(gpsTime))

    # return elements from Time Series
    return fp.data.data, fc.data.data
  else: # we'll have to work out the value at each point in the time series
    fp = np.zeros(len(gpsTime))
    fc = np.zeros(len(gpsTime))

    for i in range(len(gpsTime)):
      gps = lal.LIGOTimeGPS( gpsTime[i] )
      gmst_rad = lal.GreenwichMeanSiderealTime(gps)

      # actual computation of antenna factors
      fp[i], fc[i] = lal.ComputeDetAMResponse(response, ra, dec, psi, gmst_rad)
  fp = fp[0]
  fc = fc[0]
  return fp, fc

