def get_filter_coefs(det,fs=4096):

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
    	print 'Detector can only be H1 or L1'

    # notch filter coefficients:
    for notchf in notchesAbsolute:
        bn, an = iir_bandstops(np.array([[notchf,1,0.1]]), fs, order=4)
        coefs.append((bn,an))

    return coefs

# and then define the filter function:
def filter_data(data_in,coefs):
    data = data_in.copy()
    for coef in coefs:
        b,a = coef
        # filtfilt applies a linear filter twice, once forward and once backwards.
        # The combined filter has linear phase.
        data = filtfilt(b, a, data)
    return data
