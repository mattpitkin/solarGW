
def iir_bandstops(fstops, fs, order=4):
    import numpy as np
    from scipy.signal import iirdesign, zpk2tf, freqz
    nyq = 0.5 * fs
    # Zeros zd, poles pd, and gain kd for the digital filter
    zd = np.array([])
    pd = np.array([])
    kd = 1

    # Notches
    for fstopData in fstops:
        fstop = fstopData[0]
        df = fstopData[1]
        df2 = fstopData[2]
        low = (fstop - df) / nyq
        high = (fstop + df) / nyq
        low2 = (fstop - df2) / nyq
        high2 = (fstop + df2) / nyq
        z, p, k = iirdesign([low,high], [low2,high2], gpass=1, gstop=6,
                            ftype='ellip', output='zpk')
        zd = np.append(zd,z)
        pd = np.append(pd,p)

    # Set gain to one at 100 Hz...better not notch there
    bPrelim,aPrelim = zpk2tf(zd, pd, 1)
    outFreq, outg0 = freqz(bPrelim, aPrelim, 100/nyq)

    # Return the numerator and denominator of the digital filter
    b,a = zpk2tf(zd,pd,k)
    return b, a

def get_filter_coefs(det,fs=4096):
    from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
    import numpy as np
    from notchfilt import iir_bandstops

    # assemble the filter b,a coefficients:
    coefs = []

    # bandpass filter parameters
    lowcut = 100
    highcut = 150
    order = 4

    # bandpass filter coefficients
    nyq = 0.5*fs
    low = lowcut / nyq
    high = highcut / nyq
    bb, ab = butter(order, [low, high], btype='band')
    coefs.append((bb,ab))

    # Frequencies of notches at known instrumental spectral line frequencies.
    # You can see these lines in the ASD above, so it is straightforward to make this list.
    if det=='L1':
    	notchesAbsolute = np.array([120.0, 139.94, 145.06, 108.992])
    elif det=='H1':
	notchesAbsolute = np.array([120.0, 139.95, 140.41, 108.992])
    else:
    	print 'Error: Detector can only be H1 or L1'
		exit()
    # notch filter coefficients:
    for notchf in notchesAbsolute:
        bn, an = iir_bandstops(np.array([[notchf,1,0.1]]), fs, order=4)
        coefs.append((bn,an))

    return coefs

# and then define the filter function:
def filter_data(data_in,coefs):
    from scipy.signal import filtfilt
    data = data_in.copy()
    for coef in coefs:
        b,a = coef
        # filtfilt applies a linear filter twice, once forward and once backwards.
        # The combined filter has linear phase.
        data = filtfilt(b, a, data)
    return data
