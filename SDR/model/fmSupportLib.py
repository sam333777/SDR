

import numpy as np
import math, cmath


def fmDemodArctan(I, Q, prev_phase = 0.0):

	# the default prev_phase phase is assumed to be zero, however
	# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

# custom function for DFT that can be used by the PSD estimate
def DFT(x):

	# number of samples
	N = len(x)

	# frequency bins
	Xf = np.zeros(N, dtype='complex')

	# iterate through all frequency bins/samples
	for m in range(N):
		for k in range(N):
			Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

	# return the vector that holds the frequency bins
	return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method used in some packages
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

	# rename the NFFT argument (notation consistent with matplotlib.psd)
	# to freq_bins (i.e., frequency bins for which we compute the spectrum)
	freq_bins = NFFT
	# frequency increment (or resolution of the frequency bins)
	df = Fs / freq_bins

	# create the frequency vector to be used on the X axis
	# for plotting the PSD on the Y axis (only positive freq)
	freq = np.arange(0, Fs / 2, df)

	# design the Hann window used to smoothen the discrete data in order
	# to reduce the spectral leakage after the Fourier transform
	hann = np.empty(freq_bins)
	for i in range(len(hann)):
		hann[i] = 0.5 * (1 - math.cos(2 * math.pi * i / (freq_bins - 1)))

	# create an empty list where the PSD for each segment is computed
	psd_list = []

	# samples should be a multiple of frequency bins, so
	# the number of segments used for estimation is an integer
	# note: for this to work you must provide an argument for the
	# number of frequency bins not greater than the number of samples!
	no_segments = int(math.floor(len(samples) / float(freq_bins)))

	# iterate through all the segments
	for k in range(no_segments):

		# apply the hann window (using pointwise multiplication)
		# before computing the Fourier transform on a segment
		windowed_samples = samples[k * freq_bins:(k + 1) * freq_bins] * hann

		# compute the Fourier transform using the built-in FFT from numpy
		Xf = np.fft.fft(windowed_samples, freq_bins)

		# since input is real, we keep only the positive half of the spectrum
		# however, we will also add the signal energy of negative frequencies
		# to have a better and more accurate PSD estimate when plotting
		Xf = Xf[0:int(freq_bins / 2)] # keep only positive freq bins
		psd_seg = (1 / (Fs * freq_bins / 2)) * (abs(Xf)**2) # compute signal power
		psd_seg = 2 * psd_seg # add the energy from the negative freq bins

		# append to the list where PSD for each segment is stored
		# in sequential order (first segment, followed by the second one, ...)
		psd_list.extend(psd_seg)

	# iterate through all the frequency bins (positive freq only)
	# from all segments and average them (one bin at a time ...)
	psd_seg = np.zeros(int(freq_bins / 2))
	for k in range(int(freq_bins / 2)):
		# iterate through all the segments
		for l in range(no_segments):
			psd_seg[k] += psd_list[k + l * int(freq_bins / 2)]
		# compute the estimate for each bin
		psd_seg[k] = psd_seg[k] / no_segments

	# translate to the decibel (dB) scale
	psd_est = np.zeros(int(freq_bins / 2))
	for k in range(int(freq_bins / 2)):
		psd_est[k] = 10 * math.log10(psd_seg[k])

	# the frequency vector and PSD estimate
	return freq, psd_est

# custom function to format the plotting of the PSD
def fmPlotPSD(ax, samples, Fs, height, title):

	# adjust grid lines as needed
	x_major_interval = (Fs / 12)
	x_minor_interval = (Fs / 12) / 4
	y_major_interval = 20
	x_epsilon = 1e-3
	# adjust x/y range as needed
	x_max = x_epsilon + Fs / 2
	x_min = 0
	y_max = 10
	y_min = y_max - 100 * height

	ax.psd(samples, NFFT=512, Fs=Fs)
	#
	# below is the custom PSD estimate, which is based on the Bartlett method
	# it is less accurate than the PSD from matplotlib, however it is sufficient
	# to help us visualize the power spectra on the acquired/filtered data
	#
	# freq, my_psd = estimatePSD(samples, NFFT=512, Fs=Fs)
	# ax.plot(freq, my_psd)
	#
	ax.set_xlim([x_min, x_max])
	ax.set_ylim([y_min, y_max])
	ax.set_xticks(np.arange(x_min, x_max, x_major_interval))
	ax.set_xticks(np.arange(x_min, x_max, x_minor_interval), minor=True)
	ax.set_yticks(np.arange(y_min, y_max, y_major_interval))
	ax.grid(which='major', alpha=0.75)
	ax.grid(which='minor', alpha=0.25)
	ax.set_xlabel('Frequency (kHz)')
	ax.set_ylabel('PSD (db/Hz)')
	ax.set_title(title)

##############################################################
# New code as part of benchmarking/testing and the project
##############################################################

# custom function to estimate PSD using the matrix approach
def matrixPSD(samples, NFFT, Fs):

	freq_bins = NFFT
	df = Fs / freq_bins
	freq = np.arange(0, Fs / 2, df)
	no_segments = int(math.floor(len(samples) / float(freq_bins)))

	# generate the DFT matrix for the given size N
	dft_matrix = np.empty((freq_bins, freq_bins), dtype='complex')
	for m in range(freq_bins):
		for k in range(freq_bins):
			dft_matrix[m, k] = cmath.exp(1j * 2 * math.pi * ((-k) * m) / freq_bins)

	# generate the Hann window for the given size N
	hann_window = np.empty(freq_bins, dtype='float')
	for i in range(freq_bins):
		hann_window[i] = 0.5 * (1 - math.cos(2 * math.pi * i / (freq_bins - 1)))

	# apply Hann window and perform matrix multiplication using nested loops
	Xf = np.zeros((no_segments, freq_bins), dtype='complex')
	for seg in range(no_segments):
		for m in range(freq_bins):
			for k in range(freq_bins):
				Xf[seg][m] += samples[seg * freq_bins + k] * hann_window[k] * dft_matrix[m][k]

	# compute power, keep only positive frequencies, average across segments, and convert to dB
	psd_est = np.zeros(int(freq_bins / 2))  # same as (freq_bins // 2)
	for m in range(freq_bins // 2):
		sum_power = 0.0
		for seg in range(no_segments):
			sum_power += (1 / ((Fs / 2) * (freq_bins / 2))) * (abs(Xf[seg][m]) ** 2)
		psd_est[m] += 10 * math.log10(sum_power / no_segments)

	return freq, psd_est

# function to unit test PSD estimation
def psdUnitTest(min=-1, max=1, Fs=1e3, size=1024, NFFT=128):

	# generate random samples for testing
	samples = np.random.uniform(low=min, high=max, size=size)

	# calculate reference PSD
	freq_ref, psd_ref = estimatePSD(samples, NFFT, Fs)

	# calculate PSD using the matrix-based function
	freq_mat, psd_mat = matrixPSD(samples, NFFT, Fs)

	# check if all the values are close within the given tolerance
	if not np.allclose(freq_ref, freq_mat, atol=1e-4):
		print("Comparison between reference frequency vectors fails")

	if not np.allclose(psd_ref, psd_mat, atol=1e-4):
		print("Comparison between reference estimate PSD and matrix PSD fails")
		print("Reference PSD:", psd_ref)
		print("Matrix PSD   :", psd_mat)
		print("Maximum difference:", np.max(np.abs(psd_ref - psd_mat)))
	else:
		print(f"Unit test for matrix PSD transform passed.")

if __name__ == "__main__":

	'''
	# this unit test (when uncommented) will confirm that
	# estimate PSD and matrix PSD are equivalent to each other
	psdUnitTest()
	'''

	# do nothing when this module is launched on its own
	pass
